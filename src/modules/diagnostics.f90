!*******************************************************************************
!< Module of diagnostics
!*******************************************************************************
module diagnostics
    use constants, only: i1, i4, i8, sp, dp
    use mhd_config_module, only: mhd_config
    use simulation_setup_module, only: ndim_field
    use mhd_data_parallel, only: nfields, ngrads
    use particle_module, only: particle_type, ptls, escaped_ptls, &
        nptl_current, nptl_escaped, nptl_escaped_max, nptl_max, &
        spherical_coord_flag, leak, leak_negp, nptl_split, &
        get_interp_paramters, get_interp_paramters_spherical, &
        pmin, pmax, COUNT_FLAG_INBOX, COUNT_FLAG_OTHERS, &
        COUNT_FLAG_ESCAPE_LX, COUNT_FLAG_ESCAPE_HX, &
        COUNT_FLAG_ESCAPE_LY, COUNT_FLAG_ESCAPE_HY, &
        COUNT_FLAG_ESCAPE_LZ, COUNT_FLAG_ESCAPE_HZ

    use mpi_module
    use hdf5
    implicit none
    private
    save
    public distributions_diagnostics, quick_check, &
        init_particle_distributions, free_particle_distributions, &
        get_pmax_global, dump_particles, &
        init_escaped_particles, free_escaped_particles, &
        reset_escaped_particles, dump_escaped_particles, &
        read_diagnostics_params

    !< Parameters for particle distributions
    real(dp) :: pmin_log, pmax_log, dp_log, dmu
    real(dp) :: pmin1, pmax1
    real(dp) :: pmin2, pmax2
    real(dp) :: pmin3, pmax3
    real(dp) :: pmin4, pmax4
    real(dp) :: dx_diag1, dy_diag1, dz_diag1
    real(dp) :: dx_diag2, dy_diag2, dz_diag2
    real(dp) :: dx_diag3, dy_diag3, dz_diag3
    real(dp) :: dx_diag4, dy_diag4, dz_diag4
    real(dp) :: pmin1_log, pmax1_log, dp1_log
    real(dp) :: pmin2_log, pmax2_log, dp2_log
    real(dp) :: pmin3_log, pmax3_log, dp3_log
    real(dp) :: pmin4_log, pmax4_log, dp4_log
    real(dp) :: dmu1, dmu2, dmu3, dmu4
    integer :: npp_global, nmu_global
    !< Only dump local distributions every few MHD output intervals
    !< When dump_interval is > # of MHD outputs, it will not dump the distribution
    !< Each local distribution can have different number of npbins and nmu
    !< Each local distribution can control their reduce factor through rx, ry, rz
    logical :: dump_local_dist1, dump_local_dist2, dump_local_dist3, dump_local_dist4
    integer :: dump_interval1, npbins1, nmu1, rx1, ry1, rz1
    integer :: dump_interval2, npbins2, nmu2, rx2, ry2, rz2
    integer :: dump_interval3, npbins3, nmu3, rx3, ry3, rz3
    integer :: dump_interval4, npbins4, nmu4, rx4, ry4, rz4
    integer :: nrx1_mhd, nry1_mhd, nrz1_mhd
    integer :: nrx2_mhd, nry2_mhd, nrz2_mhd
    integer :: nrx3_mhd, nry3_mhd, nrz3_mhd
    integer :: nrx4_mhd, nry4_mhd, nrz4_mhd
    integer :: nrx1, nry1, nrz1
    integer :: nrx2, nry2, nrz2
    integer :: nrx3, nry3, nrz3
    integer :: nrx4, nry4, nrz4

    !< 2D (1D for Parker transport) distributions of p and mu
    real(dp), allocatable, dimension(:, :) :: fglobal, fglobal_sum
    real(dp), allocatable, dimension(:) :: pbins_edges_global
    real(dp), allocatable, dimension(:) :: mubins_edges_global

    !< Local distributions (mu, p, x, y, z)
    real(dp), allocatable, dimension(:) :: pbins1_edges, mubins1_edges
    real(dp), allocatable, dimension(:) :: pbins2_edges, mubins2_edges
    real(dp), allocatable, dimension(:) :: pbins3_edges, mubins3_edges
    real(dp), allocatable, dimension(:) :: pbins4_edges, mubins4_edges
    real(dp), allocatable, dimension(:, :, :, :, :) :: flocal1, flocal1_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: flocal2, flocal2_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: flocal3, flocal3_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: flocal4, flocal4_sum

    !< Distributions of p and mu for the escaped particles
    !< They use the same pbins_edges_global and mubins_edges_global as above
    real(dp), allocatable, dimension(:, :, :) :: fescaped, fescaped_sum

    !< Local distributions for the escaped particles
    !< They use the same pbins and mubins as above
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped1_x, fescaped1_x_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped1_y, fescaped1_y_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped1_z, fescaped1_z_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped2_x, fescaped2_x_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped2_y, fescaped2_y_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped2_z, fescaped2_z_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped3_x, fescaped3_x_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped3_y, fescaped3_y_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped3_z, fescaped3_z_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped4_x, fescaped4_x_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped4_y, fescaped4_y_sum
    real(dp), allocatable, dimension(:, :, :, :, :) :: fescaped4_z, fescaped4_z_sum


    interface write_ptl_element
        module procedure &
            write_integer1_element, write_integer4_element, &
            write_integer8_element, write_double_element
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
        integer(i8) :: i
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
    !< Args:
    !<  local_dist: whether to dump local particle distributions
    !<  dump_escaped_dist: whether escaped particle distributions
    !---------------------------------------------------------------------------
    subroutine init_particle_distributions(local_dist, dump_escaped_dist)
        implicit none
        logical, intent(in) :: local_dist, dump_escaped_dist
        integer :: i
        allocate(fglobal(nmu_global, npp_global))
        if (mpi_cross_rank == master) then
            allocate(fglobal_sum(nmu_global, npp_global))
        endif

        !< escaped particles
        if (dump_escaped_dist) then
            allocate(fescaped(nmu_global, npp_global, ndim_field*2))
            if (mpi_cross_rank == master) then
                allocate(fescaped_sum(nmu_global, npp_global, ndim_field*2))
            endif
        endif

        !< Momentum bins and mu bins. These are the edges of the bins.
        allocate(pbins_edges_global(npp_global+1))
        allocate(mubins_edges_global(nmu_global+1))
        pbins_edges_global = 0.0_dp
        pmin_log = log10(pmin)
        pmax_log = log10(pmax)
        dp_log = (pmax_log - pmin_log) / npp_global
        do i = 1, npp_global+1
            pbins_edges_global(i) = 10**(pmin_log + (i-1) * dp_log)
        enddo
        mubins_edges_global = 0.0_dp
        dmu = 2.0 / nmu_global
        do i = 1, nmu_global+1
            mubins_edges_global(i) = -1.0 + (i - 1) * dmu
        enddo

        if (local_dist) then
            call init_local_particle_distributions
            if (dump_escaped_dist) then
                call init_local_escaped_distributions
            endif
        endif

        call clean_particle_distributions(local_dist)
        if (dump_escaped_dist) then
            call clean_escaped_distributions(local_dist)
        endif

        if (mpi_rank == master) then
            write(*, "(A)") "Finished Initializing particle distributions."
        endif
    end subroutine init_particle_distributions

    !---------------------------------------------------------------------------
    !< Initialize local particle distributions
    !---------------------------------------------------------------------------
    subroutine init_local_particle_distributions
        implicit none
        integer :: i
        if (dump_local_dist1) then
            allocate(flocal1(nmu1, npbins1, nrx1, nry1, nrz1))
        endif
        if (dump_local_dist2) then
            allocate(flocal2(nmu2, npbins2, nrx2, nry2, nrz2))
        endif
        if (dump_local_dist3) then
            allocate(flocal3(nmu3, npbins3, nrx3, nry3, nrz3))
        endif
        if (dump_local_dist4) then
            allocate(flocal4(nmu4, npbins4, nrx4, nry4, nrz4))
        endif

        if (mpi_cross_rank == master) then
            if (dump_local_dist1) then
                allocate(flocal1_sum(nmu1, npbins1, nrx1, nry1, nrz1))
            endif
            if (dump_local_dist2) then
                allocate(flocal2_sum(nmu2, npbins2, nrx2, nry2, nrz2))
            endif
            if (dump_local_dist3) then
                allocate(flocal3_sum(nmu3, npbins3, nrx3, nry3, nrz3))
            endif
            if (dump_local_dist4) then
                allocate(flocal4_sum(nmu4, npbins4, nrx4, nry4, nrz4))
            endif
        endif

        !< Parameters for local distributions
        if (dump_local_dist1) then
            nrx1_mhd = (mhd_config%nx + rx1 - 1) / rx1
            nry1_mhd = (mhd_config%ny + ry1 - 1) / ry1
            nrz1_mhd = (mhd_config%nz + rz1 - 1) / rz1
            dx_diag1 = mhd_config%lx / nrx1_mhd
            dy_diag1 = mhd_config%ly / nry1_mhd
            dz_diag1 = mhd_config%lz / nrz1_mhd
            pmin1_log = log10(pmin1)
            pmax1_log = log10(pmax1)
            dp1_log = (pmax1_log - pmin1_log) / npbins1
            dmu1 = 2.0 / nmu1
            allocate(pbins1_edges(npbins1+1))
            allocate(mubins1_edges(nmu1+1))
            pbins1_edges = 0.0_dp
            do i = 1, npbins1+1
                pbins1_edges(i) = 10**(pmin1_log + (i-1) * dp1_log)
            enddo
            mubins1_edges = 0.0_dp
            do i = 1, nmu1+1
                mubins1_edges(i) = -1.0 + (i - 1) * dmu1
            enddo
        endif

        if (dump_local_dist2) then
            nrx2_mhd = (mhd_config%nx + rx2 - 1) / rx2
            nry2_mhd = (mhd_config%ny + ry2 - 1) / ry2
            nrz2_mhd = (mhd_config%nz + rz2 - 1) / rz2
            dx_diag2 = mhd_config%lx / nrx2_mhd
            dy_diag2 = mhd_config%ly / nry2_mhd
            dz_diag2 = mhd_config%lz / nrz2_mhd
            pmin2_log = log10(pmin2)
            pmax2_log = log10(pmax2)
            dp2_log = (pmax2_log - pmin2_log) / npbins2
            dmu2 = 2.0 / nmu2
            allocate(pbins2_edges(npbins2+1))
            allocate(mubins2_edges(nmu2+1))
            pbins2_edges = 0.0_dp
            do i = 1, npbins2+1
                pbins2_edges(i) = 10**(pmin2_log + (i-1) * dp2_log)
            enddo
            mubins2_edges = 0.0_dp
            do i = 1, nmu2+1
                mubins2_edges(i) = -1.0 + (i - 1) * dmu2
            enddo
        endif

        if (dump_local_dist3) then
            nrx3_mhd = (mhd_config%nx + rx3 - 1) / rx3
            nry3_mhd = (mhd_config%ny + ry3 - 1) / ry3
            nrz3_mhd = (mhd_config%nz + rz3 - 1) / rz3
            dx_diag3 = mhd_config%lx / nrx3_mhd
            dy_diag3 = mhd_config%ly / nry3_mhd
            dz_diag3 = mhd_config%lz / nrz3_mhd
            pmin3_log = log10(pmin3)
            pmax3_log = log10(pmax3)
            dp3_log = (pmax3_log - pmin3_log) / npbins3
            dmu3 = 2.0 / nmu3
            allocate(pbins3_edges(npbins3+1))
            allocate(mubins3_edges(nmu3+1))
            pbins3_edges = 0.0_dp
            do i = 1, npbins3+1
                pbins3_edges(i) = 10**(pmin3_log + (i-1) * dp3_log)
            enddo
            mubins3_edges = 0.0_dp
            do i = 1, nmu3+1
                mubins3_edges(i) = -1.0 + (i - 1) * dmu3
            enddo
        endif

        if (dump_local_dist4) then
            nrx4_mhd = (mhd_config%nx + rx4 - 1) / rx4
            nry4_mhd = (mhd_config%ny + ry4 - 1) / ry4
            nrz4_mhd = (mhd_config%nz + rz4 - 1) / rz4
            dx_diag4 = mhd_config%lx / nrx4_mhd
            dy_diag4 = mhd_config%ly / nry4_mhd
            dz_diag4 = mhd_config%lz / nrz4_mhd
            pmin4_log = log10(pmin4)
            pmax4_log = log10(pmax4)
            dp4_log = (pmax4_log - pmin4_log) / npbins4
            dmu4 = 2.0 / nmu4
            allocate(pbins4_edges(npbins4+1))
            allocate(mubins4_edges(nmu4+1))
            pbins4_edges = 0.0_dp
            do i = 1, npbins4+1
                pbins4_edges(i) = 10**(pmin4_log + (i-1) * dp4_log)
            enddo
            mubins4_edges = 0.0_dp
            do i = 1, nmu4+1
                mubins4_edges(i) = -1.0 + (i - 1) * dmu4
            enddo
        endif
    end subroutine init_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Initialize local distributions for escaped particles
    !---------------------------------------------------------------------------
    subroutine init_local_escaped_distributions
        implicit none
        integer :: i
        if (dump_local_dist1) then
            allocate(fescaped1_x(nmu1, npbins1, nry1, nrz1, 2))
            if (ndim_field > 1) then
                allocate(fescaped1_y(nmu1, npbins1, nrx1, nrz1, 2))
            endif
            if (ndim_field > 2) then
                allocate(fescaped1_z(nmu1, npbins1, nrx1, nry1, 2))
            endif
        endif
        if (dump_local_dist2) then
            allocate(fescaped2_x(nmu2, npbins2, nry2, nrz2, 2))
            if (ndim_field > 1) then
                allocate(fescaped2_y(nmu2, npbins2, nrx2, nrz2, 2))
            endif
            if (ndim_field > 2) then
                allocate(fescaped2_z(nmu2, npbins2, nrx2, nry2, 2))
            endif
        endif
        if (dump_local_dist3) then
            allocate(fescaped3_x(nmu3, npbins3, nry3, nrz3, 2))
            if (ndim_field > 1) then
                allocate(fescaped3_y(nmu3, npbins3, nrx3, nrz3, 2))
            endif
            if (ndim_field > 2) then
                allocate(fescaped3_z(nmu3, npbins3, nrx3, nry3, 2))
            endif
        endif
        if (dump_local_dist4) then
            allocate(fescaped4_x(nmu4, npbins4, nry4, nrz4, 2))
            if (ndim_field > 1) then
                allocate(fescaped4_y(nmu4, npbins4, nrx4, nrz4, 2))
            endif
            if (ndim_field > 2) then
                allocate(fescaped4_z(nmu4, npbins4, nrx4, nry4, 2))
            endif
        endif

        if (mpi_cross_rank == master) then
            if (dump_local_dist1) then
                allocate(fescaped1_x_sum(nmu1, npbins1, nry1, nrz1, 2))
                if (ndim_field > 1) then
                    allocate(fescaped1_y_sum(nmu1, npbins1, nrx1, nrz1, 2))
                endif
                if (ndim_field > 2) then
                    allocate(fescaped1_z_sum(nmu1, npbins1, nrx1, nry1, 2))
                endif
            endif
            if (dump_local_dist2) then
                allocate(fescaped2_x_sum(nmu2, npbins2, nry2, nrz2, 2))
                if (ndim_field > 1) then
                    allocate(fescaped2_y_sum(nmu2, npbins2, nrx2, nrz2, 2))
                endif
                if (ndim_field > 2) then
                    allocate(fescaped2_z_sum(nmu2, npbins2, nrx2, nry2, 2))
                endif
            endif
            if (dump_local_dist3) then
                allocate(fescaped3_x_sum(nmu3, npbins3, nry3, nrz3, 2))
                if (ndim_field > 1) then
                    allocate(fescaped3_y_sum(nmu3, npbins3, nrx3, nrz3, 2))
                endif
                if (ndim_field > 2) then
                    allocate(fescaped3_z_sum(nmu3, npbins3, nrx3, nry3, 2))
                endif
            endif
            if (dump_local_dist4) then
                allocate(fescaped4_x_sum(nmu4, npbins4, nry4, nrz4, 2))
                if (ndim_field > 1) then
                    allocate(fescaped4_y_sum(nmu4, npbins4, nrx4, nrz4, 2))
                endif
                if (ndim_field > 2) then
                    allocate(fescaped4_z_sum(nmu4, npbins4, nrx4, nry4, 2))
                endif
            endif
        endif
    end subroutine init_local_escaped_distributions

    !---------------------------------------------------------------------------
    !< Set particle distributions to be zero
    !< Args:
    !<  local_dist: whether to dump local particle distributions
    !---------------------------------------------------------------------------
    subroutine clean_particle_distributions(local_dist)
        implicit none
        logical, intent(in) :: local_dist
        fglobal = 0.0_dp
        if (mpi_cross_rank == master) then
            fglobal_sum = 0.0_dp
        endif
        if (local_dist) then
            call clean_local_particle_distribution
        endif
    end subroutine clean_particle_distributions

    !---------------------------------------------------------------------------
    !< Set local particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_local_particle_distribution
        implicit none
        if (dump_local_dist1) then
            flocal1 = 0.0_dp
            if (mpi_cross_rank == master) then
                flocal1_sum = 0.0_dp
            endif
        endif
        if (dump_local_dist2) then
            flocal2 = 0.0_dp
            if (mpi_cross_rank == master) then
                flocal2_sum = 0.0_dp
            endif
        endif
        if (dump_local_dist3) then
            flocal3 = 0.0_dp
            if (mpi_cross_rank == master) then
                flocal3_sum = 0.0_dp
            endif
        endif
        if (dump_local_dist4) then
            flocal4 = 0.0_dp
            if (mpi_cross_rank == master) then
                flocal4_sum = 0.0_dp
            endif
        endif
    end subroutine clean_local_particle_distribution

    !---------------------------------------------------------------------------
    !< Set the distributions of escaped particles to be zero
    !< Args:
    !<  local_dist: whether to dump local particle distributions
    !---------------------------------------------------------------------------
    subroutine clean_escaped_distributions(local_dist)
        implicit none
        logical, intent(in) :: local_dist
        fescaped = 0.0_dp
        if (mpi_cross_rank == master) then
            fescaped_sum = 0.0_dp
        endif
        if (local_dist) then
            call clean_local_escaped_distribution
        endif
    end subroutine clean_escaped_distributions

    !---------------------------------------------------------------------------
    !< Set local distributions for escaped particles to be zero
    !---------------------------------------------------------------------------
    subroutine clean_local_escaped_distribution
        implicit none
        if (dump_local_dist1) then
            fescaped1_x = 0.0_dp
            if (ndim_field > 1) then
                fescaped1_y = 0.0_dp
            endif
            if (ndim_field > 2) then
                fescaped1_z = 0.0_dp
            endif
            if (mpi_cross_rank == master) then
                fescaped1_x_sum = 0.0_dp
                if (ndim_field > 1) then
                    fescaped1_y_sum = 0.0_dp
                endif
                if (ndim_field > 2) then
                    fescaped1_z_sum = 0.0_dp
                endif
            endif
        endif

        if (dump_local_dist2) then
            fescaped2_x = 0.0_dp
            if (ndim_field > 1) then
                fescaped2_y = 0.0_dp
            endif
            if (ndim_field > 2) then
                fescaped2_z = 0.0_dp
            endif
            if (mpi_cross_rank == master) then
                fescaped2_x_sum = 0.0_dp
                if (ndim_field > 1) then
                    fescaped2_y_sum = 0.0_dp
                endif
                if (ndim_field > 2) then
                    fescaped2_z_sum = 0.0_dp
                endif
            endif
        endif

        if (dump_local_dist3) then
            fescaped3_x = 0.0_dp
            if (ndim_field > 1) then
                fescaped3_y = 0.0_dp
            endif
            if (ndim_field > 2) then
                fescaped3_z = 0.0_dp
            endif
            if (mpi_cross_rank == master) then
                fescaped3_x_sum = 0.0_dp
                if (ndim_field > 1) then
                    fescaped3_y_sum = 0.0_dp
                endif
                if (ndim_field > 2) then
                    fescaped3_z_sum = 0.0_dp
                endif
            endif
        endif

        if (dump_local_dist4) then
            fescaped4_x = 0.0_dp
            if (ndim_field > 1) then
                fescaped4_y = 0.0_dp
            endif
            if (ndim_field > 2) then
                fescaped4_z = 0.0_dp
            endif
            if (mpi_cross_rank == master) then
                fescaped4_x_sum = 0.0_dp
                if (ndim_field > 1) then
                    fescaped4_y_sum = 0.0_dp
                endif
                if (ndim_field > 2) then
                    fescaped4_z_sum = 0.0_dp
                endif
            endif
        endif
    end subroutine clean_local_escaped_distribution

    !---------------------------------------------------------------------------
    !< Free particle distributions
    !< Args:
    !<  local_dist: whether to dump local particle distributions
    !<  dump_escaped_dist: whether escaped particle distributions
    !---------------------------------------------------------------------------
    subroutine free_particle_distributions(local_dist, dump_escaped_dist)
        implicit none
        logical, intent(in) :: local_dist, dump_escaped_dist
        deallocate(fglobal)
        deallocate(pbins_edges_global)
        deallocate(mubins_edges_global)
        if (mpi_cross_rank == master) then
            deallocate(fglobal_sum)
        endif

        if (dump_escaped_dist) then
            deallocate(fescaped)
            if (mpi_cross_rank == master) then
                deallocate(fescaped_sum)
            endif
        endif

        if (local_dist) then
            call free_local_particle_distributions
            if (dump_escaped_dist) then
                call free_local_escaped_distributions
            endif
        endif
    end subroutine free_particle_distributions

    !---------------------------------------------------------------------------
    !< Free local particle distributions
    !---------------------------------------------------------------------------
    subroutine free_local_particle_distributions
        implicit none
        if (dump_local_dist1) then
            deallocate(flocal1, pbins1_edges, mubins1_edges)
            if (mpi_cross_rank == master) then
                deallocate(flocal1_sum)
            endif
        endif

        if (dump_local_dist2) then
            deallocate(flocal2, pbins2_edges, mubins2_edges)
            if (mpi_cross_rank == master) then
                deallocate(flocal2_sum)
            endif
        endif

        if (dump_local_dist3) then
            deallocate(flocal3, pbins3_edges, mubins3_edges)
            if (mpi_cross_rank == master) then
                deallocate(flocal3_sum)
            endif
        endif

        if (dump_local_dist4) then
            deallocate(flocal4, pbins4_edges, mubins4_edges)
            if (mpi_cross_rank == master) then
                deallocate(flocal4_sum)
            endif
        endif
    end subroutine free_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Free the distributions for escaped particles
    !---------------------------------------------------------------------------
    subroutine free_local_escaped_distributions
        implicit none
        if (dump_local_dist1) then
            deallocate(fescaped1_x)
            if (ndim_field > 1) then
                deallocate(fescaped1_y)
            endif
            if (ndim_field > 2) then
                deallocate(fescaped1_z)
            endif
            if (mpi_cross_rank == master) then
                deallocate(fescaped1_x_sum)
                if (ndim_field > 1) then
                    deallocate(fescaped1_y_sum)
                endif
                if (ndim_field > 2) then
                    deallocate(fescaped1_z_sum)
                endif
            endif
        endif

        if (dump_local_dist2) then
            deallocate(fescaped2_x)
            if (ndim_field > 1) then
                deallocate(fescaped2_y)
            endif
            if (ndim_field > 2) then
                deallocate(fescaped2_z)
            endif
            if (mpi_cross_rank == master) then
                deallocate(fescaped2_x_sum)
                if (ndim_field > 1) then
                    deallocate(fescaped2_y_sum)
                endif
                if (ndim_field > 2) then
                    deallocate(fescaped2_z_sum)
                endif
            endif
        endif

        if (dump_local_dist3) then
            deallocate(fescaped3_x)
            if (ndim_field > 1) then
                deallocate(fescaped3_y)
            endif
            if (ndim_field > 2) then
                deallocate(fescaped3_z)
            endif
            if (mpi_cross_rank == master) then
                deallocate(fescaped3_x_sum)
                if (ndim_field > 1) then
                    deallocate(fescaped3_y_sum)
                endif
                if (ndim_field > 2) then
                    deallocate(fescaped3_z_sum)
                endif
            endif
        endif

        if (dump_local_dist4) then
            deallocate(fescaped4_x)
            if (ndim_field > 1) then
                deallocate(fescaped4_y)
            endif
            if (ndim_field > 2) then
                deallocate(fescaped4_z)
            endif
            if (mpi_cross_rank == master) then
                deallocate(fescaped4_x_sum)
                if (ndim_field > 1) then
                    deallocate(fescaped4_y_sum)
                endif
                if (ndim_field > 2) then
                    deallocate(fescaped4_z_sum)
                endif
            endif
        endif
    end subroutine free_local_escaped_distributions

    !---------------------------------------------------------------------------
    !< Accumulate particle distributions
    !< Args:
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine calc_particle_distributions(local_dist)
        use simulation_setup_module, only: fconfig
        implicit none
        logical, intent(in) :: local_dist
        integer(i8) :: iptl
        integer :: ip, imu
        integer :: ix1, iy1, iz1, ip1, imu1
        integer :: ix2, iy2, iz2, ip2, imu2
        integer :: ix3, iy3, iz3, ip3, imu3
        integer :: ix4, iy4, iz4, ip4, imu4
        real(dp) :: xmin, ymin, zmin
        real(dp) :: x, y, z, p, mu, weight
        logical :: condx, condy, condz, condp, condmu
        type(particle_type) :: ptl

        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin

        do iptl = 1, nptl_current
            ptl = ptls(iptl)
            x = ptl%x
            y = ptl%y
            z = ptl%z
            p = ptl%p
            mu = ptl%mu
            weight = ptl%weight
            if (spherical_coord_flag) then
                ! For spherical coordinates, we solve for F=f*p^2*r^2 for 1D and
                ! f*p^2*r^2*sin(theta) for 2D and 3D
                weight = weight / (ptl%x)**2
                if (ndim_field >= 2) then
                    weight = weight / sin(ptl%y)
                endif
            endif

            ! Global distributions
            if (p > pmin .and. p <= pmax .and. &
                mu >= -1.0d0 .and. mu <= 1.0d0) then
                ip = floor((log10(p)-pmin_log) / dp_log) + 1
                imu = floor((mu + 1.0d0) / dmu) + 1
                fglobal(imu,ip) = fglobal(imu,ip) + ptl%weight
            endif

            if (local_dist) then
                if (dump_local_dist1) then
                    ! Local distributions 1
                    ix1 = floor((x - xmin)/dx_diag1) + 1
                    iy1 = floor((y - ymin)/dy_diag1) + 1
                    iz1 = floor((z - zmin)/dz_diag1) + 1
                    ip1 = floor((log10(p)-pmin1_log) / dp1_log) + 1
                    imu1 = floor((mu + 1.0d0) / dmu1) + 1
                    condx = ix1 >= 1 .and. ix1 <= nrx1
                    condy = iy1 >= 1 .and. iy1 <= nry1
                    condz = iz1 >= 1 .and. iz1 <= nrz1
                    condp = ip1 > 0 .and. ip1 < npbins1
                    condmu = imu1 >=1 .and. imu1 <= nmu1
                    if (condx .and. condy .and. condz .and. condp .and. condmu) then
                        flocal1(imu1, ip1, ix1, iy1, iz1) = &
                            flocal1(imu1, ip1, ix1, iy1, iz1) + weight
                    endif
                endif

                if (dump_local_dist2) then
                    ! Local distributions 2
                    ix2 = floor((x - xmin)/dx_diag2) + 1
                    iy2 = floor((y - ymin)/dy_diag2) + 1
                    iz2 = floor((z - zmin)/dz_diag2) + 1
                    ip2 = floor((log10(p)-pmin2_log) / dp2_log) + 1
                    imu2 = floor((mu + 1.0d0) / dmu2) + 1
                    condx = ix2 >= 1 .and. ix2 <= nrx2
                    condy = iy2 >= 1 .and. iy2 <= nry2
                    condz = iz2 >= 1 .and. iz2 <= nrz2
                    condp = ip2 > 0 .and. ip2 < npbins2
                    condmu = imu2 >=1 .and. imu2 <= nmu2
                    if (condx .and. condy .and. condz .and. condp .and. condmu) then
                        flocal2(imu2, ip2, ix2, iy2, iz2) = &
                            flocal2(imu2, ip2, ix2, iy2, iz2) + weight
                    endif
                endif

                if (dump_local_dist3) then
                    ! Local distributions 3
                    ix3 = floor((x - xmin)/dx_diag3) + 1
                    iy3 = floor((y - ymin)/dy_diag3) + 1
                    iz3 = floor((z - zmin)/dz_diag3) + 1
                    ip3 = floor((log10(p)-pmin3_log) / dp3_log) + 1
                    imu3 = floor((mu + 1.0d0) / dmu3) + 1
                    condx = ix3 >= 1 .and. ix3 <= nrx3
                    condy = iy3 >= 1 .and. iy3 <= nry3
                    condz = iz3 >= 1 .and. iz3 <= nrz3
                    condp = ip3 > 0 .and. ip3 < npbins3
                    condmu = imu3 >=1 .and. imu3 <= nmu3
                    if (condx .and. condy .and. condz .and. condp .and. condmu) then
                        flocal3(imu3, ip3, ix3, iy3, iz3) = &
                            flocal3(imu3, ip3, ix3, iy3, iz3) + weight
                    endif
                endif

                if (dump_local_dist4) then
                    ! Local distributions 4
                    ix4 = floor((x - xmin)/dx_diag4) + 1
                    iy4 = floor((y - ymin)/dy_diag4) + 1
                    iz4 = floor((z - zmin)/dz_diag4) + 1
                    ip4 = floor((log10(p)-pmin4_log) / dp4_log) + 1
                    imu4 = floor((mu + 1.0d0) / dmu4) + 1
                    condx = ix4 >= 1 .and. ix4 <= nrx4
                    condy = iy4 >= 1 .and. iy4 <= nry4
                    condz = iz4 >= 1 .and. iz4 <= nrz4
                    condp = ip4 > 0 .and. ip4 < npbins4
                    condmu = imu4 >=1 .and. imu4 <= nmu4
                    if (condx .and. condy .and. condz .and. condp .and. condmu) then
                        flocal4(imu4, ip4, ix4, iy4, iz4) = &
                            flocal4(imu4, ip4, ix4, iy4, iz4) + weight
                    endif
                endif
            endif
        enddo ! Loop over particles

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fglobal, fglobal_sum, nmu_global*npp_global, &
            MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
        if (local_dist) then
            if (dump_local_dist1) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(flocal1, flocal1_sum, nmu1*npbins1*nrx1*nry1*nrz1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
            endif
            if (dump_local_dist2) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(flocal2, flocal2_sum, nmu2*npbins2*nrx2*nry2*nrz2, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
            endif
            if (dump_local_dist3) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(flocal3, flocal3_sum, nmu3*npbins3*nrx3*nry3*nrz3, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
            endif
            if (dump_local_dist4) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(flocal4, flocal4_sum, nmu4*npbins4*nrx4*nry4*nrz4, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
            endif
        endif
    end subroutine calc_particle_distributions

    !---------------------------------------------------------------------------
    !< Accumulate the distributions for escaped particles
    !< Args:
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine calc_escaped_distributions(local_dist)
        use simulation_setup_module, only: fconfig
        implicit none
        logical, intent(in) :: local_dist
        integer(i8) :: iptl
        integer :: ip, imu
        integer :: ix1, iy1, iz1, ip1, imu1
        integer :: ix2, iy2, iz2, ip2, imu2
        integer :: ix3, iy3, iz3, ip3, imu3
        integer :: ix4, iy4, iz4, ip4, imu4
        integer :: ibc
        real(dp) :: xmin, ymin, zmin
        real(dp) :: x, y, z, p, mu, weight
        logical :: condx, condy, condz, condp, condmu
        type(particle_type) :: ptl

        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin

        do iptl = 1, nptl_escaped
            ptl = escaped_ptls(iptl)
            x = ptl%x
            y = ptl%y
            z = ptl%z
            p = ptl%p
            mu = ptl%mu
            weight = ptl%weight
            if (spherical_coord_flag) then
                ! For spherical coordinates, we solve for F=f*p^2*r^2 for 1D and
                ! f*p^2*r^2*sin(theta) for 2D and 3D
                weight = weight / (ptl%x)**2
                if (ndim_field >= 2) then
                    weight = weight / sin(ptl%y)
                endif
            endif

            ! Global distributions
            if (p > pmin .and. p <= pmax .and. &
                mu >= -1.0d0 .and. mu <= 1.0d0) then
                ip = floor((log10(p)-pmin_log) / dp_log) + 1
                imu = floor((mu + 1.0d0) / dmu) + 1
                ibc = abs(ptl%count_flag)
                fescaped(imu,ip,ibc) = fescaped(imu,ip,ibc) + ptl%weight
            endif

            if (local_dist) then
                if (dump_local_dist1) then
                    ! Local distributions 1
                    ix1 = floor((x - xmin)/dx_diag1) + 1
                    iy1 = floor((y - ymin)/dy_diag1) + 1
                    iz1 = floor((z - zmin)/dz_diag1) + 1
                    ip1 = floor((log10(p)-pmin1_log) / dp1_log) + 1
                    imu1 = floor((mu + 1.0d0) / dmu1) + 1
                    condx = ix1 >= 1 .and. ix1 <= nrx1
                    condy = iy1 >= 1 .and. iy1 <= nry1
                    condz = iz1 >= 1 .and. iz1 <= nrz1
                    condp = ip1 > 0 .and. ip1 < npbins1
                    condmu = imu1 >=1 .and. imu1 <= nmu1
                    if (condp .and. condmu) then
                        if (ptl%count_flag == COUNT_FLAG_ESCAPE_LX) then
                            if (condy .and. condz) then
                                fescaped1_x(imu1, ip1, iy1, iz1, 1) = &
                                    fescaped1_x(imu1, ip1, iy1, iz1, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HX) then
                            if (condy .and. condz) then
                                fescaped1_x(imu1, ip1, iy1, iz1, 2) = &
                                    fescaped1_x(imu1, ip1, iy1, iz1, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LY) then
                            if (condx .and. condz) then
                                fescaped1_y(imu1, ip1, ix1, iz1, 1) = &
                                    fescaped1_y(imu1, ip1, ix1, iz1, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HY) then
                            if (condx .and. condz) then
                                fescaped1_z(imu1, ip1, ix1, iz1, 2) = &
                                    fescaped1_z(imu1, ip1, ix1, iz1, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LZ) then
                            if (condx .and. condy) then
                                fescaped1_y(imu1, ip1, ix1, iy1, 1) = &
                                    fescaped1_y(imu1, ip1, ix1, iy1, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HZ) then
                            if (condx .and. condy) then
                                fescaped1_z(imu1, ip1, ix1, iy1, 2) = &
                                    fescaped1_z(imu1, ip1, ix1, iy1, 2) + weight
                            endif
                        endif
                    endif
                endif

                if (dump_local_dist2) then
                    ! Local distributions 2
                    ix2 = floor((x - xmin)/dx_diag2) + 1
                    iy2 = floor((y - ymin)/dy_diag2) + 1
                    iz2 = floor((z - zmin)/dz_diag2) + 1
                    ip2 = floor((log10(p)-pmin2_log) / dp2_log) + 1
                    imu2 = floor((mu + 1.0d0) / dmu2) + 1
                    condx = ix2 >= 1 .and. ix2 <= nrx2
                    condy = iy2 >= 1 .and. iy2 <= nry2
                    condz = iz2 >= 1 .and. iz2 <= nrz2
                    condp = ip2 > 0 .and. ip2 < npbins2
                    condmu = imu2 >=1 .and. imu2 <= nmu2
                    if (condp .and. condmu) then
                        if (ptl%count_flag == COUNT_FLAG_ESCAPE_LX) then
                            if (condy .and. condz) then
                                fescaped2_x(imu2, ip2, iy2, iz2, 1) = &
                                    fescaped2_x(imu2, ip2, iy2, iz2, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HX) then
                            if (condy .and. condz) then
                                fescaped2_x(imu2, ip2, iy2, iz2, 2) = &
                                    fescaped2_x(imu2, ip2, iy2, iz2, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LY) then
                            if (condx .and. condz) then
                                fescaped2_y(imu2, ip2, ix2, iz2, 1) = &
                                    fescaped2_y(imu2, ip2, ix2, iz2, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HY) then
                            if (condx .and. condz) then
                                fescaped2_z(imu2, ip2, ix2, iz2, 2) = &
                                    fescaped2_z(imu2, ip2, ix2, iz2, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LZ) then
                            if (condx .and. condy) then
                                fescaped2_y(imu2, ip2, ix2, iy2, 1) = &
                                    fescaped2_y(imu2, ip2, ix2, iy2, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HZ) then
                            if (condx .and. condy) then
                                fescaped2_z(imu2, ip2, ix2, iy2, 2) = &
                                    fescaped2_z(imu2, ip2, ix2, iy2, 2) + weight
                            endif
                        endif
                    endif
                endif

                if (dump_local_dist3) then
                    ! Local distributions 3
                    ix3 = floor((x - xmin)/dx_diag3) + 1
                    iy3 = floor((y - ymin)/dy_diag3) + 1
                    iz3 = floor((z - zmin)/dz_diag3) + 1
                    ip3 = floor((log10(p)-pmin3_log) / dp3_log) + 1
                    imu3 = floor((mu + 1.0d0) / dmu3) + 1
                    condx = ix3 >= 1 .and. ix3 <= nrx3
                    condy = iy3 >= 1 .and. iy3 <= nry3
                    condz = iz3 >= 1 .and. iz3 <= nrz3
                    condp = ip3 > 0 .and. ip3 < npbins3
                    condmu = imu3 >=1 .and. imu3 <= nmu3
                    if (condp .and. condmu) then
                        if (ptl%count_flag == COUNT_FLAG_ESCAPE_LX) then
                            if (condy .and. condz) then
                                fescaped3_x(imu3, ip3, iy3, iz3, 1) = &
                                    fescaped3_x(imu3, ip3, iy3, iz3, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HX) then
                            if (condy .and. condz) then
                                fescaped3_x(imu3, ip3, iy3, iz3, 2) = &
                                    fescaped3_x(imu3, ip3, iy3, iz3, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LY) then
                            if (condx .and. condz) then
                                fescaped3_y(imu3, ip3, ix3, iz3, 1) = &
                                    fescaped3_y(imu3, ip3, ix3, iz3, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HY) then
                            if (condx .and. condz) then
                                fescaped3_z(imu3, ip3, ix3, iz3, 2) = &
                                    fescaped3_z(imu3, ip3, ix3, iz3, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LZ) then
                            if (condx .and. condy) then
                                fescaped3_y(imu3, ip3, ix3, iy3, 1) = &
                                    fescaped3_y(imu3, ip3, ix3, iy3, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HZ) then
                            if (condx .and. condy) then
                                fescaped3_z(imu3, ip3, ix3, iy3, 2) = &
                                    fescaped3_z(imu3, ip3, ix3, iy3, 2) + weight
                            endif
                        endif
                    endif
                endif

                if (dump_local_dist4) then
                    ! Local distributions 4
                    ix4 = floor((x - xmin)/dx_diag4) + 1
                    iy4 = floor((y - ymin)/dy_diag4) + 1
                    iz4 = floor((z - zmin)/dz_diag4) + 1
                    ip4 = floor((log10(p)-pmin4_log) / dp4_log) + 1
                    imu4 = floor((mu + 1.0d0) / dmu4) + 1
                    condx = ix4 >= 1 .and. ix4 <= nrx4
                    condy = iy4 >= 1 .and. iy4 <= nry4
                    condz = iz4 >= 1 .and. iz4 <= nrz4
                    condp = ip4 > 0 .and. ip4 < npbins4
                    condmu = imu4 >=1 .and. imu4 <= nmu4
                    if (condp .and. condmu) then
                        if (ptl%count_flag == COUNT_FLAG_ESCAPE_LX) then
                            if (condy .and. condz) then
                                fescaped4_x(imu4, ip4, iy4, iz4, 1) = &
                                    fescaped4_x(imu4, ip4, iy4, iz4, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HX) then
                            if (condy .and. condz) then
                                fescaped4_x(imu4, ip4, iy4, iz4, 2) = &
                                    fescaped4_x(imu4, ip4, iy4, iz4, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LY) then
                            if (condx .and. condz) then
                                fescaped4_y(imu4, ip4, ix4, iz4, 1) = &
                                    fescaped4_y(imu4, ip4, ix4, iz4, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HY) then
                            if (condx .and. condz) then
                                fescaped4_z(imu4, ip4, ix4, iz4, 2) = &
                                    fescaped4_z(imu4, ip4, ix4, iz4, 2) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_LZ) then
                            if (condx .and. condy) then
                                fescaped4_y(imu4, ip4, ix4, iy4, 1) = &
                                    fescaped4_y(imu4, ip4, ix4, iy4, 1) + weight
                            endif
                        else if (ptl%count_flag == COUNT_FLAG_ESCAPE_HZ) then
                            if (condx .and. condy) then
                                fescaped4_z(imu4, ip4, ix4, iy4, 2) = &
                                    fescaped4_z(imu4, ip4, ix4, iy4, 2) + weight
                            endif
                        endif
                    endif
                endif
            endif
        enddo ! Loop over particles

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fescaped, fescaped_sum, nmu_global*npp_global*ndim_field*2, &
            MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
        if (local_dist) then
            if (dump_local_dist1) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(fescaped1_x, fescaped1_x_sum, nmu1*npbins1*nry1*nrz1*2, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                if (ndim_field > 1) then
                    call MPI_REDUCE(fescaped1_y, fescaped1_y_sum, nmu1*npbins1*nrx1*nrz1*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
                if (ndim_field > 2) then
                    call MPI_REDUCE(fescaped1_z, fescaped1_z_sum, nmu1*npbins1*nrx1*nry1*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
            endif
            if (dump_local_dist2) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(fescaped2_x, fescaped2_x_sum, nmu2*npbins2*nry2*nrz2*2, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                if (ndim_field > 1) then
                    call MPI_REDUCE(fescaped2_y, fescaped2_y_sum, nmu2*npbins2*nrx2*nrz2*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
                if (ndim_field > 2) then
                    call MPI_REDUCE(fescaped2_z, fescaped2_z_sum, nmu2*npbins2*nrx2*nry2*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
            endif
            if (dump_local_dist3) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(fescaped3_x, fescaped3_x_sum, nmu3*npbins3*nry3*nrz3*2, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                if (ndim_field > 1) then
                    call MPI_REDUCE(fescaped3_y, fescaped3_y_sum, nmu3*npbins3*nrx3*nrz3*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
                if (ndim_field > 2) then
                    call MPI_REDUCE(fescaped3_z, fescaped3_z_sum, nmu3*npbins3*nrx3*nry3*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
            endif
            if (dump_local_dist4) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_REDUCE(fescaped4_x, fescaped4_x_sum, nmu4*npbins4*nry4*nrz4*2, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                if (ndim_field > 1) then
                    call MPI_REDUCE(fescaped4_y, fescaped4_y_sum, nmu4*npbins4*nrx4*nrz4*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
                if (ndim_field > 2) then
                    call MPI_REDUCE(fescaped4_z, fescaped4_z_sum, nmu4*npbins4*nrx4*nry4*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
                endif
            endif
        endif
    end subroutine calc_escaped_distributions

    !---------------------------------------------------------------------------
    !< Save the momentum and mu bins for the global distributions
    !< Args:
    !<  file_id: file ID for the HDF5 file
    !<  nmu: number of mu bins
    !<  npp: number of momentum bins
    !<  dname_mubins: dataset name for mu bins
    !<  dname_pbins: dataset name for momentum bins
    !<  mubins_edges: the edges of the mu bins
    !<  pbins_edges: the edges of the momentum bins
    !---------------------------------------------------------------------------
    subroutine save_mu_pbins(file_id, nmu, npp, dname_mubins, dname_pbins, &
            mubins_edges, pbins_edges)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer, intent(in) :: nmu, npp
        character(*), intent(in) :: dname_mubins, dname_pbins
        real(dp), dimension(:), intent(in) :: mubins_edges, pbins_edges
        integer(hid_t) :: dset_id, filespace
        integer(hsize_t), dimension(1) :: dcount_1d, doffset_1d, dset_dims_1d
        integer :: error

        ! Mu bins edges
        dcount_1d(1) = nmu + 1
        doffset_1d(1) = 0
        dset_dims_1d(1) = nmu + 1
        call h5screate_simple_f(1, dset_dims_1d, filespace, error)
        call h5dcreate_f(file_id, trim(dname_mubins), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount_1d, doffset_1d, dset_dims_1d, mubins_edges)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)

        ! Momentum bins edges
        dcount_1d(1) = npp + 1
        doffset_1d(1) = 0
        dset_dims_1d(1) = npp + 1
        call h5screate_simple_f(1, dset_dims_1d, filespace, error)
        call h5dcreate_f(file_id, trim(dname_pbins), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount_1d, doffset_1d, dset_dims_1d, pbins_edges)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine save_mu_pbins

    !---------------------------------------------------------------------------
    !< Save global particle distributions using HDF5 format
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_global_distributions(iframe, file_path)
        use mpi_io_module, only: set_mpi_datatype_double, set_mpi_info, &
            fileinfo, open_data_mpi_io, write_data_mpi_io
        use hdf5_io, only: create_file_h5, close_file_h5, write_data_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        character(len=4) :: ctime
        character(len=128) :: fname
        integer(hid_t) :: file_id, dset_id, filespace
        integer(hsize_t), dimension(2) :: dcount_2d, doffset_2d, dset_dims_2d
        integer :: error

        call h5open_f(error)

        write (ctime,'(i4.4)') iframe
        if (mpi_rank == master) then
            fname = trim(file_path)//'fdists_'//ctime//'.h5'
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .false., mpi_sub_comm)
            call save_mu_pbins(file_id, nmu_global, npp_global, &
                "mubins_edges_global", "pbins_edges_global", &
                mubins_edges_global, pbins_edges_global)
            ! Distribution data
            dcount_2d(1) = nmu_global
            dcount_2d(2) = npp_global
            doffset_2d(1) = 0
            doffset_2d(2) = 0
            dset_dims_2d(1) = nmu_global
            dset_dims_2d(2) = npp_global
            call h5screate_simple_f(2, dset_dims_2d, filespace, error)
            call h5dcreate_f(file_id, "fglobal", H5T_NATIVE_DOUBLE, &
                filespace, dset_id, error)
            call write_data_h5(dset_id, dcount_2d, doffset_2d, dset_dims_2d, fglobal_sum)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(filespace, error)
            call close_file_h5(file_id)
        endif
        call h5close_f(error)
    end subroutine save_global_distributions

    !---------------------------------------------------------------------------
    !< Save the distributions for escaped particles using HDF5 format
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_escaped_distributions(iframe, file_path)
        use mpi_io_module, only: set_mpi_datatype_double, set_mpi_info, &
            fileinfo, open_data_mpi_io, write_data_mpi_io
        use hdf5_io, only: create_file_h5, close_file_h5, write_data_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        character(len=4) :: ctime
        character(len=128) :: fname
        integer(hid_t) :: file_id, dset_id, filespace
        integer(hsize_t), dimension(3) :: dcount_3d, doffset_3d, dset_dims_3d
        integer :: error

        call h5open_f(error)

        write (ctime,'(i4.4)') iframe
        if (mpi_rank == master) then
            fname = trim(file_path)//'escaped_dists_'//ctime//'.h5'
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .false., mpi_sub_comm)
            call save_mu_pbins(file_id, nmu_global, npp_global, &
                "mubins_edges", "pbins_edges", &
                mubins_edges_global, pbins_edges_global)
            ! Distribution data
            dcount_3d(1) = nmu_global
            dcount_3d(2) = npp_global
            dcount_3d(3) = ndim_field*2
            doffset_3d(1) = 0
            doffset_3d(2) = 0
            doffset_3d(3) = 0
            dset_dims_3d(1) = nmu_global
            dset_dims_3d(2) = npp_global
            dset_dims_3d(3) = ndim_field*2
            call h5screate_simple_f(3, dset_dims_3d, filespace, error)
            call h5dcreate_f(file_id, "fescaped", H5T_NATIVE_DOUBLE, &
                filespace, dset_id, error)
            call write_data_h5(dset_id, dcount_3d, doffset_3d, dset_dims_3d, fescaped_sum)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(filespace, error)
            call close_file_h5(file_id)
        endif
        call h5close_f(error)
    end subroutine save_escaped_distributions

    !---------------------------------------------------------------------------
    !< Save the 5D particle distributions data
    !< Args:
    !<  fname: filename for the HDF5 file
    !<  nmu: number of mu bins
    !<  npbins: number of momentum bins
    !<  nrx, nry, nrz: number of cells (reduced) along each dimension locally
    !<  nrx_mhd: number of cells (reduced) along x globally
    !<  nry_mhd: number of cells (reduced) along y globally
    !<  nrz_mhd: number of cells (reduced) along z globally
    !<  dset_name: the dataset name
    !<  flocal: the data for the local distributions
    !---------------------------------------------------------------------------
    subroutine save_5d_local_dist_data(fname, nmu, npbins, nrx, nry, nrz, &
            nrx_mhd, nry_mhd, nrz_mhd, dset_name, flocal)
        use hdf5_io, only: open_file_h5, close_file_h5, write_data_h5
        use simulation_setup_module, only: mpi_ix, mpi_iy, mpi_iz
        implicit none
        character(*), intent(in) :: fname, dset_name
        integer, intent(in) :: nmu, npbins, nrx, nry, nrz
        integer, intent(in) :: nrx_mhd, nry_mhd, nrz_mhd
        real(dp), dimension(:, :, :, :, :), intent(in) :: flocal
        integer(hid_t) :: file_id, dset_id, filespace
        integer(hsize_t), dimension(5) :: dcount_5d, doffset_5d, dset_dims_5d
        integer :: error

        call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .true., mpi_sub_comm)
        dcount_5d = (/ nmu, npbins, nrx, nry, nrz /)
        doffset_5d = (/ 0, 0, nrx * mpi_ix, nry * mpi_iy, nrz * mpi_iz /)
        dset_dims_5d = (/ nmu, npbins, nrx_mhd, nry_mhd, nrz_mhd /)
        call h5screate_simple_f(5, dset_dims_5d, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount_5d, doffset_5d, dset_dims_5d, &
            flocal, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
        call close_file_h5(file_id)
    end subroutine save_5d_local_dist_data

    !---------------------------------------------------------------------------
    !< Save local particle distributions using HDF5 format
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_local_distributions(iframe, file_path)
        use hdf5_io, only: open_file_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        character(len=4) :: ctime
        character(len=128) :: fname
        integer(hid_t) :: file_id
        integer :: error

        call h5open_f(error)

        write (ctime,'(i4.4)') iframe
        if (mpi_cross_rank == master) then
            fname = trim(file_path)//'fdists_'//ctime//'.h5'
            ! Local distributions 1
            if (dump_local_dist1) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu1, npbins1, &
                        "mubins1_edges", "pbins1_edges", &
                        mubins1_edges, pbins1_edges)
                    call close_file_h5(file_id)
                endif
                call save_5d_local_dist_data(fname, nmu1, npbins1, nrx1, nry1, nrz1, &
                    nrx1_mhd, nry1_mhd, nrz1_mhd, "flocal1", flocal1_sum)
            endif

            ! Local distributions 2
            if (dump_local_dist2) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu2, npbins2, &
                        "mubins2_edges", "pbins2_edges", &
                        mubins2_edges, pbins2_edges)
                    call close_file_h5(file_id)
                endif
                call save_5d_local_dist_data(fname, nmu2, npbins2, nrx2, nry2, nrz2, &
                    nrx2_mhd, nry2_mhd, nrz2_mhd, "flocal2", flocal2_sum)
            endif

            ! Local distributions 3
            if (dump_local_dist3) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu3, npbins3, &
                        "mubins3_edges", "pbins3_edges", &
                        mubins3_edges, pbins3_edges)
                    call close_file_h5(file_id)
                endif
                call save_5d_local_dist_data(fname, nmu3, npbins3, nrx3, nry3, nrz3, &
                    nrx3_mhd, nry3_mhd, nrz3_mhd, "flocal3", flocal3_sum)
            endif

            ! Local distributions 4
            if (dump_local_dist4) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu4, npbins4, &
                        "mubins4_edges", "pbins4_edges", &
                        mubins4_edges, pbins4_edges)
                    call close_file_h5(file_id)
                endif
                call save_5d_local_dist_data(fname, nmu4, npbins4, nrx4, nry4, nrz4, &
                    nrx4_mhd, nry4_mhd, nrz4_mhd, "flocal4", flocal4_sum)
            endif
        endif
        call h5close_f(error)
    end subroutine save_local_distributions

    !---------------------------------------------------------------------------
    !< Save the 5D distributions for escaped particles data
    !< Args:
    !<  fname: filename for the HDF5 file
    !<  nmu: number of mu bins
    !<  npbins: number of momentum bins
    !<  nr1, nr2: number of cells (reduced) along each dimension locally
    !<  nr1_mhd: number of cells (reduced) along one dimension globally
    !<  nr2_mhd: number of cells (reduced) along the other dimension globally
    !<  mpi_i1: MPI rank along one dimension
    !<  mpi_i2: MPI rank along the other dimension
    !<  dset_name: the dataset name
    !<  flocal: the data for the local distributions
    !---------------------------------------------------------------------------
    subroutine save_5d_local_escaped_dist_data(fname, nmu, npbins, nr1, nr2, &
            nr1_mhd, nr2_mhd, mpi_i1, mpi_i2, dset_name, flocal)
        use hdf5_io, only: open_file_h5, close_file_h5, write_data_h5
        implicit none
        character(*), intent(in) :: fname, dset_name
        integer, intent(in) :: nmu, npbins, nr1, nr2
        integer, intent(in) :: nr1_mhd, nr2_mhd, mpi_i1, mpi_i2
        real(dp), dimension(:, :, :, :, :), intent(in) :: flocal
        integer(hid_t) :: file_id, dset_id, filespace
        integer(hsize_t), dimension(5) :: dcount_5d, doffset_5d, dset_dims_5d
        integer :: error

        call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .true., mpi_sub_comm)
        dcount_5d = (/ nmu, npbins, nr1, nr2, 2 /)
        doffset_5d = (/ 0, 0, nr1 * mpi_i1, nr2 * mpi_i2, 0 /)
        dset_dims_5d = (/ nmu, npbins, nr1_mhd, nr2_mhd, 2 /)
        call h5screate_simple_f(5, dset_dims_5d, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount_5d, doffset_5d, dset_dims_5d, &
            flocal, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
        call close_file_h5(file_id)
    end subroutine save_5d_local_escaped_dist_data

    !---------------------------------------------------------------------------
    !< Save local distributions for escaped particles using HDF5 format
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_local_escaped_distributions(iframe, file_path)
        use mpi_io_module, only: set_mpi_datatype_double, set_mpi_info, &
            fileinfo, open_data_mpi_io, write_data_mpi_io
        use hdf5_io, only: open_file_h5, close_file_h5, write_data_h5
        use simulation_setup_module, only: mpi_ix, mpi_iy, mpi_iz
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        character(len=4) :: ctime
        character(len=128) :: fname
        integer(hid_t) :: file_id
        integer :: error

        call h5open_f(error)

        write (ctime,'(i4.4)') iframe
        if (mpi_cross_rank == master) then
            fname = trim(file_path)//'escaped_dists_'//ctime//'.h5'
            ! Local distributions 1
            if (dump_local_dist1) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu1, npbins1, &
                        "mubins1_edges", "pbins1_edges", &
                        mubins1_edges, pbins1_edges)
                    call close_file_h5(file_id)
                endif
                call save_5d_local_escaped_dist_data(fname, nmu1, npbins1, nry1, nrz1, &
                        nry1_mhd, nrz1_mhd, mpi_iy, mpi_iz, "fescaped1_x", fescaped1_x_sum)
                if (ndim_field > 1) then
                    call save_5d_local_escaped_dist_data(fname, nmu1, npbins1, nrx1, nrz1, &
                            nrx1_mhd, nrz1_mhd, mpi_ix, mpi_iz, "fescaped1_y", fescaped1_y_sum)
                endif
                if (ndim_field > 2) then
                    call save_5d_local_escaped_dist_data(fname, nmu1, npbins1, nrx1, nry1, &
                            nrx1_mhd, nry1_mhd, mpi_ix, mpi_iy, "fescaped1_z", fescaped1_z_sum)
                endif
            endif

            ! Local distributions 2
            if (dump_local_dist2) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu2, npbins2, &
                        "mubins2_edges", "pbins2_edges", &
                        mubins2_edges, pbins2_edges)
                    call close_file_h5(file_id)
                endif

                call save_5d_local_escaped_dist_data(fname, nmu2, npbins2, nry2, nrz2, &
                        nry2_mhd, nrz2_mhd, mpi_iy, mpi_iz, "fescaped2_x", fescaped2_x_sum)
                if (ndim_field > 1) then
                    call save_5d_local_escaped_dist_data(fname, nmu2, npbins2, nrx2, nrz2, &
                            nrx2_mhd, nrz2_mhd, mpi_ix, mpi_iz, "fescaped2_y", fescaped2_y_sum)
                endif
                if (ndim_field > 2) then
                    call save_5d_local_escaped_dist_data(fname, nmu2, npbins2, nrx2, nry2, &
                            nrx2_mhd, nry2_mhd, mpi_ix, mpi_iy, "fescaped2_z", fescaped2_z_sum)
                endif
            endif

            ! Local distributions 3
            if (dump_local_dist3) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu3, npbins3, &
                        "mubins3_edges", "pbins3_edges", &
                        mubins3_edges, pbins3_edges)
                    call close_file_h5(file_id)
                endif

                call save_5d_local_escaped_dist_data(fname, nmu3, npbins3, nry3, nrz3, &
                        nry3_mhd, nrz3_mhd, mpi_iy, mpi_iz, "fescaped3_x", fescaped3_x_sum)
                if (ndim_field > 1) then
                    call save_5d_local_escaped_dist_data(fname, nmu3, npbins3, nrx3, nrz3, &
                            nrx3_mhd, nrz3_mhd, mpi_ix, mpi_iz, "fescaped3_y", fescaped3_y_sum)
                endif
                if (ndim_field > 2) then
                    call save_5d_local_escaped_dist_data(fname, nmu3, npbins3, nrx3, nry3, &
                            nrx3_mhd, nry3_mhd, mpi_ix, mpi_iy, "fescaped3_z", fescaped3_z_sum)
                endif
            endif

            ! Local distributions 4
            if (dump_local_dist4) then
                if (mpi_sub_rank == master) then
                    call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .false.)
                    call save_mu_pbins(file_id, nmu4, npbins4, &
                        "mubins4_edges", "pbins4_edges", &
                        mubins4_edges, pbins4_edges)
                    call close_file_h5(file_id)
                endif

                call save_5d_local_escaped_dist_data(fname, nmu4, npbins4, nry4, nrz4, &
                        nry4_mhd, nrz4_mhd, mpi_iy, mpi_iz, "fescaped4_x", fescaped4_x_sum)
                if (ndim_field > 1) then
                    call save_5d_local_escaped_dist_data(fname, nmu4, npbins4, nrx4, nrz4, &
                            nrx4_mhd, nrz4_mhd, mpi_ix, mpi_iz, "fescaped4_y", fescaped4_y_sum)
                endif
                if (ndim_field > 2) then
                    call save_5d_local_escaped_dist_data(fname, nmu4, npbins4, nrx4, nry4, &
                            nrx4_mhd, nry4_mhd, mpi_ix, mpi_iy, "fescaped4_z", fescaped4_z_sum)
                endif
            endif
        endif
        call h5close_f(error)
    end subroutine save_local_escaped_distributions

    !---------------------------------------------------------------------------
    !< Diagnostics of the particle distributions
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !<  local_dist: whether to accumulate local particle distribution
    !<  dump_escaped_dist: whether escaped particle distributions
    !---------------------------------------------------------------------------
    subroutine distributions_diagnostics(iframe, file_path, local_dist, dump_escaped_dist)
        use constants, only: sp, dp
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: local_dist, dump_escaped_dist
        character(*), intent(in) :: file_path
        integer :: fh, pos1
        character(len=4) :: ctime
        character(len=128) :: fname
        logical :: dir_e

        call calc_particle_distributions(local_dist)
        call save_global_distributions(iframe, file_path)

        if (dump_escaped_dist) then
            call calc_escaped_distributions(local_dist)
            call save_escaped_distributions(iframe, file_path)
        endif

        if (local_dist) then
            call save_local_distributions(iframe, file_path)
            if (dump_escaped_dist) then
                call save_local_escaped_distributions(iframe, file_path)
            endif
        endif

        call clean_particle_distributions(local_dist)
        if (dump_escaped_dist) then
            call clean_escaped_distributions(local_dist)
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
        real(dp) :: pmax_local, pmax_global
        integer(i8) :: i
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
    !< Write one char-length integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_integer1_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i1), dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_STD_I8LE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer1_element

    !---------------------------------------------------------------------------
    !< Write one default-length integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_integer4_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i4), dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_INTEGER, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer4_element

    !---------------------------------------------------------------------------
    !< Write one long-length integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_integer8_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i8), dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_STD_I64LE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer8_element

    !---------------------------------------------------------------------------
    !< dump all the particles
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, open_file_h5, close_file_h5, write_data_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        integer(i8) :: nptl_local, nptl_global, nptl_offset
        integer(i8), allocatable, dimension(:) :: nptls_local
        character(len=256) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id, dset_id, filespace
        integer :: error
        logical :: dir_e

        call h5open_f(error)

        ! Write the local number of particles
        write (ctime,'(i4.4)') iframe
        fname = trim(file_path)//'particles_'//ctime//'.h5'
        allocate(nptls_local(mpi_size))
        call MPI_GATHER(nptl_current, 1, MPI_INTEGER8, nptls_local(mpi_rank+1), &
            1, MPI_INTEGER8, master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .false., MPI_COMM_WORLD)
            dcount(1) = mpi_size
            doffset(1) = 0
            dset_dims(1) = mpi_size
            call h5screate_simple_f(1, dset_dims, filespace, error)
            call h5dcreate_f(file_id, "nptl_local", H5T_STD_I64LE, &
                filespace, dset_id, error)
            call write_data_h5(dset_id, dcount, doffset, dset_dims, &
                nptls_local, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(filespace, error)
            call close_file_h5(file_id)
        endif
        deallocate(nptls_local)

        ! Write the particle data
        call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .true., MPI_COMM_WORLD)
        nptl_local = nptl_current
        call MPI_ALLREDUCE(nptl_local, nptl_global, 1, MPI_INTEGER8, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_SCAN(nptl_local, nptl_offset, 1, MPI_INTEGER8, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        nptl_offset = nptl_offset - nptl_local

        dcount(1) = nptl_local
        doffset(1) = nptl_offset
        dset_dims(1) = nptl_global

        call write_ptl_element(file_id, dcount, doffset, dset_dims, "x", ptls%x)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "y", ptls%y)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "z", ptls%z)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "p", ptls%p)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "v", ptls%v)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "mu", ptls%mu)
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
            "parent", ptls%parent)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "nsteps_tracking", ptls%nsteps_tracking)

        call close_file_h5(file_id)
        call h5close_f(error)
    end subroutine dump_particles

    !---------------------------------------------------------------------------
    !< dump escaped particles
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_escaped_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, open_file_h5, close_file_h5, write_data_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        integer(i8) :: nptl_local, nptl_global, nptl_offset
        integer(i8), allocatable, dimension(:) :: nptls_local
        character(len=128) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id, dset_id, filespace
        integer :: error
        logical :: dir_e

        nptl_local = nptl_escaped
        call MPI_ALLREDUCE(nptl_local, nptl_global, 1, MPI_INTEGER8, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_SCAN(nptl_local, nptl_offset, 1, MPI_INTEGER8, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        nptl_offset = nptl_offset - nptl_local

        if (nptl_global > 0) then
            CALL h5open_f(error)

            ! Write the local number of particles
            write (ctime,'(i4.4)') iframe
            fname = trim(file_path)//'escaped_particles_'//ctime//'.h5'
            allocate(nptls_local(mpi_size))
            call MPI_GATHER(nptl_current, 1, MPI_INTEGER8, nptls_local(mpi_rank+1), &
                1, MPI_INTEGER8, master, MPI_COMM_WORLD, ierr)
            if (mpi_rank == master) then
                call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .false., MPI_COMM_WORLD)
                dcount(1) = mpi_size
                doffset(1) = 0
                dset_dims(1) = mpi_size
                call h5screate_simple_f(1, dset_dims, filespace, error)
                call h5dcreate_f(file_id, "nptl_local", H5T_STD_I64LE, &
                    filespace, dset_id, error)
                call write_data_h5(dset_id, dcount, doffset, dset_dims, &
                    nptls_local, .false., .false.)
                call h5dclose_f(dset_id, error)
                call h5sclose_f(filespace, error)
                call close_file_h5(file_id)
            endif
            deallocate(nptls_local)

            dcount(1) = nptl_local
            doffset(1) = nptl_offset
            dset_dims(1) = nptl_global

            call open_file_h5(fname, H5F_ACC_RDWR_F, file_id, .true., MPI_COMM_WORLD)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "x", escaped_ptls%x)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "y", escaped_ptls%y)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "z", escaped_ptls%z)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "p", escaped_ptls%p)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "v", escaped_ptls%v)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "mu", escaped_ptls%mu)
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
                "parent", escaped_ptls%parent)
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
    !< Check the configuration for the local spatial diagnostics
    !< Args:
    !---------------------------------------------------------------------------
    subroutine check_local_dist_configuration(nx, ny, nz, rx, ry, rz)
        use simulation_setup_module, only: fconfig
        implicit none
        integer, intent(in) :: nx, ny, nz, rx, ry, rz
        if (ndim_field == 1) then
            if (nx * rx /= fconfig%nx) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'rx' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        else if (ndim_field == 2) then
            if (nx * rx /= fconfig%nx .or. &
                ny * ry /= fconfig%ny) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'ry' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        else
            if (nx * rx /= fconfig%nx .or. &
                ny * ry /= fconfig%ny .or. &
                nz * rz /= fconfig%nz) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'rz' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        endif
    end subroutine check_local_dist_configuration

    !---------------------------------------------------------------------------
    !< Read parameters for diagnostics
    !---------------------------------------------------------------------------
    subroutine echo_diagnostics_params(dump_interval, pmin, pmax, nx, ny, nz, &
            npbins, nmu, rx, ry, rz)
        implicit none
        real(dp), intent(in) :: pmin, pmax
        integer, intent(in) :: dump_interval, nx, ny, nz, npbins, nmu, rx, ry, rz
        if (mpi_rank == master) then
            write(*, "(A)") ""
            write(*, "(A)") " Parameters for the local distribution diagnostics "
            write(*, "(A,I0,A)") "  Dump distributions every ", &
                dump_interval, " MHD frames"
            write(*, "(A, E13.6E2, E13.6E2)") "  Min and Max particle momentum: ", &
                pmin, pmax
            write(*, "(A,I0,A,I0,A,I0)") "  Dimensions: ", nx, " ", ny, " ", nz
            write(*, "(A,I0)") "  Number of momentum bins: ", npbins
            write(*, "(A,I0)") "  Number of mu bins: ", nmu
            write(*, "(A,I0)") "  Reduce factor along x: ", rx
            write(*, "(A,I0)") "  Reduce factor along y: ", ry
            write(*, "(A,I0)") "  Reduce factor along z: ", rz
        endif
    end subroutine echo_diagnostics_params

    !---------------------------------------------------------------------------
    !< Read parameters for diagnostics
    !< Args:
    !<  conf_file: configuration file name
    !<  focused_transport: whether to the Focused Transport equation
    !<  nframes: number of MHD frames
    !---------------------------------------------------------------------------
    subroutine read_diagnostics_params(conf_file, focused_transport, nframes)
        use read_config, only: get_variable
        use simulation_setup_module, only: fconfig
        implicit none
        character(*), intent(in) :: conf_file
        logical, intent(in) :: focused_transport
        integer, intent(in) :: nframes
        integer, parameter :: ndiag_params = 34
        real(dp), dimension(ndiag_params) :: diag_params
        integer :: fh

        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/'//trim(conf_file), status='old')
            diag_params(1) = get_variable(fh, 'npp_global', '=')
            diag_params(2) = get_variable(fh, 'nmu_global', '=')

            diag_params(3) = get_variable(fh, 'dump_interval1', '=')
            diag_params(4) = get_variable(fh, 'pmin1', '=')
            diag_params(5) = get_variable(fh, 'pmax1', '=')
            diag_params(6) = get_variable(fh, 'npbins1', '=')
            diag_params(7) = get_variable(fh, 'nmu1', '=')
            diag_params(8) = get_variable(fh, 'rx1', '=')
            diag_params(9) = get_variable(fh, 'ry1', '=')
            diag_params(10) = get_variable(fh, 'rz1', '=')

            diag_params(11) = get_variable(fh, 'dump_interval2', '=')
            diag_params(12) = get_variable(fh, 'pmin2', '=')
            diag_params(13) = get_variable(fh, 'pmax2', '=')
            diag_params(14) = get_variable(fh, 'npbins2', '=')
            diag_params(15) = get_variable(fh, 'nmu2', '=')
            diag_params(16) = get_variable(fh, 'rx2', '=')
            diag_params(17) = get_variable(fh, 'ry2', '=')
            diag_params(18) = get_variable(fh, 'rz2', '=')

            diag_params(19) = get_variable(fh, 'dump_interval3', '=')
            diag_params(20) = get_variable(fh, 'pmin3', '=')
            diag_params(21) = get_variable(fh, 'pmax3', '=')
            diag_params(22) = get_variable(fh, 'npbins3', '=')
            diag_params(23) = get_variable(fh, 'nmu3', '=')
            diag_params(24) = get_variable(fh, 'rx3', '=')
            diag_params(25) = get_variable(fh, 'ry3', '=')
            diag_params(26) = get_variable(fh, 'rz3', '=')

            diag_params(27) = get_variable(fh, 'dump_interval4', '=')
            diag_params(28) = get_variable(fh, 'pmin4', '=')
            diag_params(29) = get_variable(fh, 'pmax4', '=')
            diag_params(30) = get_variable(fh, 'npbins4', '=')
            diag_params(31) = get_variable(fh, 'nmu4', '=')
            diag_params(32) = get_variable(fh, 'rx4', '=')
            diag_params(33) = get_variable(fh, 'ry4', '=')
            diag_params(34) = get_variable(fh, 'rz4', '=')
            close(fh)
        endif
        call MPI_BCAST(diag_params, ndiag_params, MPI_DOUBLE_PRECISION, &
            master, MPI_COMM_WORLD, ierr)
        npp_global = int(diag_params(1))
        nmu_global = int(diag_params(2))
        if (focused_transport) then
            nmu_global = int(diag_params(2))
        else
            nmu_global = 1 ! For Parker's transport
        endif
        if (mpi_rank == master) then
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A,I0)") " # of momentum bins for global distribution: ", npp_global
            write(*, "(A,I0)") " # of mu bins for global distribution: ", nmu_global
        endif

        ! 1st set of local distribution diagnostics parameters
        dump_interval1 = int(diag_params(3))
        pmin1 = diag_params(4)
        pmax1 = diag_params(5)
        npbins1 = int(diag_params(6))
        if (focused_transport) then
            nmu1 = int(diag_params(7))
        else
            nmu1 = 1 ! For Parker's transport
        endif
        if (dump_interval1 < nframes) then ! Otherwise, the set of parameters will be used
            rx1 = int(diag_params(8))
            ry1 = int(diag_params(9))
            rz1 = int(diag_params(10))
            nrx1 = (fconfig%nx + rx1 - 1) / rx1
            nry1 = (fconfig%ny + ry1 - 1) / ry1
            nrz1 = (fconfig%nz + rz1 - 1) / rz1
            call echo_diagnostics_params(dump_interval1, pmin1, pmax1, nrx1, nry1, nrz1, &
                                         npbins1, nmu1, rx1, ry1, rz1)
            call check_local_dist_configuration(nrx1, nry1, nrz1, rx1, ry1, rz1)
            dump_local_dist1 = .true.
        else
            dump_local_dist1 = .false.
        endif

        ! 2nd set of local distribution diagnostics parameters
        dump_interval2 = int(diag_params(11))
        pmin2 = diag_params(12)
        pmax2 = diag_params(13)
        npbins2 = int(diag_params(14))
        if (focused_transport) then
            nmu2 = int(diag_params(15))
        else
            nmu2 = 1 ! For Parker's transport
        endif
        if (dump_interval2 < nframes) then ! Otherwise, the set of parameters will be used
            rx2 = int(diag_params(16))
            ry2 = int(diag_params(17))
            rz2 = int(diag_params(18))
            nrx2 = (fconfig%nx + rx2 - 1) / rx2
            nry2 = (fconfig%ny + ry2 - 1) / ry2
            nrz2 = (fconfig%nz + rz2 - 1) / rz2
            call echo_diagnostics_params(dump_interval2, pmin2, pmax2, nrx2, nry2, nrz2, &
                                         npbins2, nmu2, rx2, ry2, rz2)
            call check_local_dist_configuration(nrx2, nry2, nrz2, rx2, ry2, rz2)
            dump_local_dist2 = .true.
        else
            dump_local_dist2 = .false.
        endif

        ! 3rd set of local distribution diagnostics parameters
        dump_interval3 = int(diag_params(19))
        pmin3 = diag_params(20)
        pmax3 = diag_params(21)
        npbins3 = int(diag_params(22))
        if (focused_transport) then
            nmu3 = int(diag_params(23))
        else
            nmu3 = 1 ! For Parker's transport
        endif
        if (dump_interval3 < nframes) then ! Otherwise, the set of parameters will be used
            rx3 = int(diag_params(24))
            ry3 = int(diag_params(25))
            rz3 = int(diag_params(26))
            nrx3 = (fconfig%nx + rx3 - 1) / rx3
            nry3 = (fconfig%ny + ry3 - 1) / ry3
            nrz3 = (fconfig%nz + rz3 - 1) / rz3
            call echo_diagnostics_params(dump_interval3, pmin3, pmax3, nrx3, nry3, nrz3, &
                                         npbins3, nmu3, rx3, ry3, rz3)
            call check_local_dist_configuration(nrx3, nry3, nrz3, rx3, ry3, rz3)
            dump_local_dist3 = .true.
        else
            dump_local_dist3 = .false.
        endif

        ! 4th set of local distribution diagnostics parameters
        dump_interval4 = int(diag_params(27))
        pmin4 = diag_params(28)
        pmax4 = diag_params(29)
        npbins4 = int(diag_params(30))
        if (focused_transport) then
            nmu4 = int(diag_params(31))
        else
            nmu4 = 1 ! For Parker's transport
        endif
        if (dump_interval4 < nframes) then ! Otherwise, the set of parameters will be used
            rx4 = int(diag_params(32))
            ry4 = int(diag_params(33))
            rz4 = int(diag_params(34))
            nrx4 = (fconfig%nx + rx4 - 1) / rx4
            nry4 = (fconfig%ny + ry4 - 1) / ry4
            nrz4 = (fconfig%nz + rz4 - 1) / rz4
            call echo_diagnostics_params(dump_interval4, pmin4, pmax4, nrx4, nry4, nrz4, &
                                         npbins4, nmu4, rx4, ry4, rz4)
            call check_local_dist_configuration(nrx4, nry4, nrz4, rx4, ry4, rz4)
            dump_local_dist4 = .true.
        else
            dump_local_dist4 = .false.
        endif

        if (mpi_rank == master) then
            print *, "---------------------------------------------------"
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
        escaped_ptls%v = 0.0
        escaped_ptls%mu = 0.0
        escaped_ptls%weight = 0.0
        escaped_ptls%t = 0.0
        escaped_ptls%dt = 0.0
        escaped_ptls%split_times = 0
        escaped_ptls%count_flag = COUNT_FLAG_OTHERS
        escaped_ptls%tag = 0
        escaped_ptls%parent = 0
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

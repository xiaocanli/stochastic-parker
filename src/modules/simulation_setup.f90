!*******************************************************************************
!< Module of simulation setup
!*******************************************************************************
module simulation_setup_module
    use mhd_config_module, only: mhd_config
    use constants, only: fp, dp
    use mpi_module
    implicit none
    private
    save

    public mpi_sizex, mpi_sizey, mpi_sizez, fconfig
    public read_simuation_mpi_topology, set_field_configuration

    integer :: mpi_sizex, mpi_sizey, mpi_sizez

    !< Dimensions of v and B in memory
    type field_configuration
        integer :: nxf, nyf, nzf
        integer :: ix_min, ix_max, iy_min, iy_max, iz_min, iz_max
    end type field_configuration

    type(field_configuration) :: fconfig

    contains

    !---------------------------------------------------------------------------
    !< Read simulation MPI topology
    !---------------------------------------------------------------------------
    subroutine read_simuation_mpi_topology
        use read_config, only: get_variable
        implicit none
        real(fp) :: temp
        integer, dimension(3) :: mpi_topology
        integer :: nx, ny, nz, fh
        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/conf.dat', status='old')
            temp = get_variable(fh, 'mpi_sizex', '=')
            mpi_topology(1) = int(temp)
            temp = get_variable(fh, 'mpi_sizey', '=')
            mpi_topology(2) = int(temp)
            temp = get_variable(fh, 'mpi_sizez', '=')
            mpi_topology(3) = int(temp)
            close(fh)
        endif

        call MPI_BCAST(mpi_topology, 3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        mpi_sizex = mpi_topology(1)
        mpi_sizey = mpi_topology(2)
        mpi_sizez = mpi_topology(3)
        nx = mhd_config%nx
        ny = mhd_config%ny
        nz = mhd_config%nz


        !< check the topology for consistency
        if (mpi_sizex*mpi_sizey*mpi_sizez /= mpi_size .or. &
            nx/mpi_sizex*mpi_sizex /= nx .or. &
            ny/mpi_sizey*mpi_sizey /= ny .or. &
            nz/mpi_sizez*mpi_sizez /= nz) then

           if (mpi_rank == master) print *, "invalid converter topology"
           call MPI_FINALIZE(ierr)
           stop

        endif
  
        if (mpi_rank == master) then
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A)") " Simulation MPI toplogy: "
            write(*, "(A,I0,A,I0,A,I0)") "x, y, z: ", mpi_sizex, " ", &
                mpi_sizey, " ", mpi_sizez
            print *, "---------------------------------------------------"
        endif
    end subroutine read_simuation_mpi_topology


    !---------------------------------------------------------------------------
    !< Set data boundaries
    !---------------------------------------------------------------------------
    subroutine set_data_boundaries(mrank, msize, n, ntot, imin, imax)
        implicit none
        integer, intent(in) :: mrank, msize, n, ntot
        integer, intent(out) :: imin, imax
        if (mrank == 0) then
            imin = 1
            if (ntot > (n + 2)) then
                imax = n + 2
            else
                imax = ntot
            endif
        else if (mrank == msize- 1) then
            imax = ntot
            if (ntot - n - 1 < 1) then
                imin = 1
            else
                imin = ntot - n - 1
            endif
        else
            imin = mrank * n - 1
            imax = (mrank + 1) * n + 2
        endif
    end subroutine set_data_boundaries

    !---------------------------------------------------------------------------
    !< Set field configuration for reading MHD data
    !< Args:
    !<  whole_data_flag: whether to load the whole dataset. 0 for no and
    !<      other numbers for yes
    !---------------------------------------------------------------------------
    subroutine set_field_configuration(whole_data_flag)
        implicit none
        integer, intent(in) :: whole_data_flag 
        integer :: ix, iy, iz, nx, ny, nz

        if (whole_data_flag == 0) then
            iz = mpi_rank / (mpi_sizex * mpi_sizey)
            iy = mod(mpi_rank, mpi_sizex * mpi_sizey) / mpi_sizex
            ix = mod(mpi_rank, mpi_sizex)

            nx = mhd_config%nx / mpi_sizex
            ny = mhd_config%ny / mpi_sizey
            nz = mhd_config%nz / mpi_sizez

            call set_data_boundaries(ix, mpi_sizex, nx, mhd_config%nx, &
                                     fconfig%ix_min, fconfig%ix_max)
            call set_data_boundaries(iy, mpi_sizey, ny, mhd_config%ny, &
                                     fconfig%iy_min, fconfig%iy_max)
            call set_data_boundaries(iz, mpi_sizez, nz, mhd_config%nz, &
                                     fconfig%iz_min, fconfig%iz_max)

            fconfig%nxf = fconfig%ix_max - fconfig%ix_min + 1
            fconfig%nyf = fconfig%iy_max - fconfig%iy_min + 1
            fconfig%nzf = fconfig%iz_max - fconfig%iz_min + 1
        else
            fconfig%ix_min = 1
            fconfig%iy_min = 1
            fconfig%iz_min = 1
            fconfig%ix_max = mhd_config%nx
            fconfig%iy_max = mhd_config%ny
            fconfig%iz_max = mhd_config%nz
            fconfig%nxf = mhd_config%nx
            fconfig%nyf = mhd_config%ny
            fconfig%nzf = mhd_config%nz
        endif
    end subroutine set_field_configuration

end module simulation_setup_module

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

    public mpi_sizex, mpi_sizey, mpi_sizez, fconfig, mpi_ix, mpi_iy, mpi_iz, &
        neighbors
    public read_simuation_mpi_topology, set_field_configuration, &
        read_particle_boundary_conditions, set_neighbors

    integer :: mpi_sizex, mpi_sizey, mpi_sizez
    integer :: mpi_ix, mpi_iy, mpi_iz
    integer :: pbcx, pbcy, pbcz
    integer, dimension(6) :: neighbors

    !< Dimensions of v and B in memory
    type field_configuration
        integer :: nxg, nyg, nzg    ! Total number of MHD grids with ghost cells
        integer :: nxf, nyf, nzf    ! Local grid sizes with ghost cells
        integer :: nx, ny, nz       ! Local grid size without ghost cells 
        integer :: ix_min, ix_max, iy_min
        integer :: iy_max, iz_min, iz_max
        real(dp) :: xmin, ymin, zmin
        real(dp) :: xmax, ymax, zmax
    end type field_configuration

    type(field_configuration) :: fconfig

    contains

    !---------------------------------------------------------------------------
    !< Read simulation MPI topology
    !< Args:
    !<  conf_file: configuration file name
    !---------------------------------------------------------------------------
    subroutine read_simuation_mpi_topology(conf_file)
        use read_config, only: get_variable
        implicit none
        character(*), intent(in) :: conf_file
        real(fp) :: temp
        integer, dimension(3) :: mpi_topology
        integer :: nx, ny, nz, fh
        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/'//trim(conf_file), status='old')
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

        mpi_iz = mpi_rank / (mpi_sizex * mpi_sizey)
        mpi_iy = mod(mpi_rank, mpi_sizex * mpi_sizey) / mpi_sizex
        mpi_ix = mod(mpi_rank, mpi_sizex)

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
            write(*, "(A,I0,A,I0,A,I0)") " x, y, z: ", mpi_sizex, " ", &
                mpi_sizey, " ", mpi_sizez
            print *, "---------------------------------------------------"
        endif
    end subroutine read_simuation_mpi_topology

    !---------------------------------------------------------------------------
    !< Read particle boundary conditions
    !< Args:
    !<  conf_file: configuration file name
    !---------------------------------------------------------------------------
    subroutine read_particle_boundary_conditions(conf_file)
        use read_config, only: get_variable
        implicit none
        character(*), intent(in) :: conf_file
        real(fp) :: temp
        integer, dimension(3) :: pbcs
        integer :: fh
        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/'//trim(conf_file), status='old')
            temp = get_variable(fh, 'pbcx', '=')
            pbcs(1) = int(temp)
            temp = get_variable(fh, 'pbcy', '=')
            pbcs(2) = int(temp)
            temp = get_variable(fh, 'pbcz', '=')
            pbcs(3) = int(temp)
            close(fh)
        endif

        call MPI_BCAST(pbcs, 3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        pbcx = pbcs(1)
        pbcy = pbcs(2)
        pbcz = pbcs(3)

        if (mpi_rank == master) then
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A)") " Particle boundary conditions: "
            if (pbcx == 0) then
                write(*, "(A)") " Particles are periodic along x."
            else
                write(*, "(A)") " Particles are open along x."
            endif
            if (pbcy == 0) then
                write(*, "(A)") " Particles are periodic along y."
            else
                write(*, "(A)") " Particles are open along y."
            endif
            if (pbcz == 0) then
                write(*, "(A)") " Particles are periodic along z."
            else
                write(*, "(A)") " Particles are open along z."
            endif
            print *, "---------------------------------------------------"
        endif
    end subroutine read_particle_boundary_conditions

    !---------------------------------------------------------------------------
    !< Set data boundaries. We assume the MHD data has been re-organized and
    !< add ghost cells on each side. Ghost cells are not included in un-resolved
    !< dimensions. For example, for 2D MHD simulation with grid nx*ny*1, the
    !< re-organized data has dimension (nx + 4) * (ny + 4) * 1
    !---------------------------------------------------------------------------
    subroutine set_data_boundaries(mrank, n, imin, imax)
        implicit none
        integer, intent(in) :: mrank, n
        integer, intent(out) :: imin, imax
        imin = mrank * n + 1
        imax = (mrank + 1) * n + 4
    end subroutine set_data_boundaries

    !---------------------------------------------------------------------------
    !< Set field configuration for reading MHD data
    !< Args:
    !<  whole_data_flag: whether to load the whole dataset. 0 for no and
    !<      other numbers for yes
    !<  ndim: the number of dimensions. 1, 2 or 3
    !---------------------------------------------------------------------------
    subroutine set_field_configuration(whole_data_flag, ndim)
        implicit none
        integer, intent(in) :: whole_data_flag, ndim 
        integer :: ix, iy, iz, nx, ny, nz

        fconfig%nxg = mhd_config%nx + 4
        if (ndim > 1) then
            fconfig%nyg = mhd_config%ny + 4
            if (ndim > 2) then
                fconfig%nzg = mhd_config%nz + 4
            else
                fconfig%nzg = mhd_config%nz
            endif
        else
            fconfig%nyg = mhd_config%ny
            fconfig%nzg = mhd_config%nz
        endif

        if (whole_data_flag == 0) then
            iz = mpi_rank / (mpi_sizex * mpi_sizey)
            iy = mod(mpi_rank, mpi_sizex * mpi_sizey) / mpi_sizex
            ix = mod(mpi_rank, mpi_sizex)

            nx = mhd_config%nx / mpi_sizex
            ny = mhd_config%ny / mpi_sizey
            nz = mhd_config%nz / mpi_sizez

            call set_data_boundaries(ix, nx, fconfig%ix_min, fconfig%ix_max)
            fconfig%xmin = ix * nx * mhd_config%dx + mhd_config%xmin
            fconfig%xmax = fconfig%xmin + nx * mhd_config%dx
            if (ndim > 1) then
                call set_data_boundaries(iy, ny, fconfig%iy_min, fconfig%iy_max)
                fconfig%ymin = iy * ny * mhd_config%dy + mhd_config%ymin
                fconfig%ymax = fconfig%ymin + ny * mhd_config%dy
                if (ndim > 2) then
                    call set_data_boundaries(iz, nz, fconfig%iz_min, fconfig%iz_max)
                    fconfig%zmin = iz * nz * mhd_config%dz + mhd_config%zmin
                    fconfig%zmax = fconfig%zmin + nz * mhd_config%dz
                else
                    fconfig%iz_min = 1
                    fconfig%iz_max = 1
                    fconfig%zmin = mhd_config%zmin
                    fconfig%zmax = mhd_config%zmax
                endif
            else
                fconfig%iy_min = 1
                fconfig%iy_max = 1
                fconfig%iz_min = 1
                fconfig%iz_max = 1
                fconfig%ymin = mhd_config%ymin
                fconfig%ymax = mhd_config%ymax
                fconfig%zmin = mhd_config%zmin
                fconfig%zmax = mhd_config%zmax
            endif

            fconfig%nxf = fconfig%ix_max - fconfig%ix_min + 1
            fconfig%nyf = fconfig%iy_max - fconfig%iy_min + 1
            fconfig%nzf = fconfig%iz_max - fconfig%iz_min + 1
            fconfig%nx = nx
            fconfig%ny = ny
            fconfig%nz = nz
        else
            fconfig%ix_min = 1
            fconfig%iy_min = 1
            fconfig%iz_min = 1
            fconfig%ix_max = mhd_config%nx + 4
            fconfig%nxf = mhd_config%nx + 4
            
            if (ndim > 1) then
                fconfig%iy_max = mhd_config%ny + 4
                fconfig%nyf = mhd_config%ny + 4
                if (ndim > 2) then
                    fconfig%iz_max = mhd_config%nz + 4
                    fconfig%nzf = mhd_config%nz + 4
                else
                    fconfig%iz_max = mhd_config%nz
                    fconfig%nzf = mhd_config%nz
                endif
            else
                fconfig%iy_max = mhd_config%ny
                fconfig%nyf = mhd_config%ny
                fconfig%iz_max = mhd_config%nz
                fconfig%nzf = mhd_config%nz
            endif

            fconfig%xmin = mhd_config%xmin
            fconfig%ymin = mhd_config%ymin
            fconfig%zmin = mhd_config%zmin

            fconfig%xmax = mhd_config%xmax
            fconfig%ymax = mhd_config%ymax
            fconfig%zmax = mhd_config%zmax

            fconfig%nx = mhd_config%nx
            fconfig%ny = mhd_config%ny
            fconfig%nz = mhd_config%nz
        endif

        if (mpi_rank == master) then
            !< echo the configuration information
            print *, "---------------------------------------------------"
            write(*, "(A)") " Field configuration in the simulation: "
            write(*, "(A,I0,A,I0,A,I0)") " nxg, nyg, nzg: ", fconfig%nxg, &
                " ", fconfig%nyg, " ", fconfig%nzg
            write(*, "(A,I0,A,I0,A,I0)") " nxf, nyf, nzf: ", fconfig%nxf, &
                " ", fconfig%nyf, " ", fconfig%nzf
            write(*, "(A,I0,A,I0,A,I0)") " nx, ny, nz: ", fconfig%nx, &
                " ", fconfig%ny, " ", fconfig%nz
            write(*, "(A,I0,A,I0,A,I0)") " ix_min, iy_min, iz_min: ", &
                fconfig%ix_min, " ", fconfig%iy_min, " ", fconfig%iz_min
            write(*, "(A,I0,A,I0,A,I0)") " ix_max, iy_max, iz_max: ", &
                fconfig%ix_max, " ", fconfig%iy_max, " ", fconfig%iz_max
            write(*, "(A,E14.7,A,E14.7,A,E14.7)") " xmin, ymin, zmin: ", &
                fconfig%xmin, " ", fconfig%ymin, " ", fconfig%zmin
            write(*, "(A,E14.7,A,E14.7,A,E14.7)") " xmax, ymax, zmax: ", &
                fconfig%xmax, " ", fconfig%ymax, " ", fconfig%zmax
            print *, "---------------------------------------------------"
        endif
    end subroutine set_field_configuration

    !---------------------------------------------------------------------------
    !< Set neighbors along one direction
    !< Args:
    !<  mrank: rank along this direction
    !<  msize: size along this direction
    !<  neighs: neighbors along this direction
    !<  pbc: particle boundary condition, 0 for periodic; 1 for open
    !---------------------------------------------------------------------------
    subroutine set_neighbor_one_direction(mrank, msize, neighs, pbc)
        implicit none
        integer, intent(in) :: mrank, msize, pbc
        integer, dimension(2), intent(out) :: neighs
        if (msize == 1) then
            if (pbc == 0) then
                neighs = mrank
            else if (pbc == 1) then
                neighs = -1
            endif
        else
            if (mrank == 0) then
                if (pbc == 0) then
                    neighs(1) = msize - 1
                else if (pbc == 1) then
                    neighs(1) = -1
                endif
                neighs(2) = 1
            else if (mrank == msize - 1) then
                neighs(1) = mrank - 1
                if (pbc == 0) then
                    neighs(2) = 0
                else if (pbc == 1) then
                    neighs(2) = -1
                endif
            else
                neighs(1) = mrank - 1
                neighs(2) = mrank + 1
            endif
        endif
    end subroutine set_neighbor_one_direction


    !---------------------------------------------------------------------------
    !< Set neighbors for each mpi_rank
    !< Args:
    !<  whole_data_flag: whether to load the whole dataset. 0 for no and
    !<      other numbers for yes
    !---------------------------------------------------------------------------
    subroutine set_neighbors(whole_data_flag)
        implicit none
        integer, intent(in) :: whole_data_flag
        integer :: ix, iy, iz, size_xy
        integer, dimension(2) :: ne(2)
        iz = mpi_rank / (mpi_sizex * mpi_sizey)
        iy = mod(mpi_rank, mpi_sizex * mpi_sizey) / mpi_sizex
        ix = mod(mpi_rank, mpi_sizex)
        size_xy = mpi_sizex * mpi_sizey
        if (whole_data_flag == 0) then
            call set_neighbor_one_direction(ix, mpi_sizex, ne, pbcx)
            if (ne(1) < 0) then
                neighbors(1) = -1
            else
                neighbors(1) = ne(1) + iy * mpi_sizex + iz * size_xy
            endif
            if (ne(2) < 0) then
                neighbors(2) = -1
            else
                neighbors(2) = ne(2) + iy * mpi_sizex + iz * size_xy
            endif

            call set_neighbor_one_direction(iy, mpi_sizey, ne, pbcy)
            if (ne(1) < 0) then
                neighbors(3) = -1
            else
                neighbors(3) = ix + ne(1) * mpi_sizex + iz * size_xy
            endif
            if (ne(2) < 0) then
                neighbors(4) = -1
            else
                neighbors(4) = ix + ne(2) * mpi_sizex + iz * size_xy
            endif

            call set_neighbor_one_direction(iz, mpi_sizez, ne, pbcz)
            if (ne(1) < 0) then
                neighbors(5) = -1
            else
                neighbors(5) = ix + iy * mpi_sizex + ne(1) * size_xy
            endif
            if (ne(2) < 0) then
                neighbors(6) = -1
            else
                neighbors(6) = ix + iy * mpi_sizex + ne(2) * size_xy
            endif
        else
            neighbors = mpi_rank
        endif
    end subroutine set_neighbors
end module simulation_setup_module

!*******************************************************************************
!< Module of MHD data for Shengtai Li's MHD code
!*******************************************************************************
module mhd_data_sli
    use constants, only: sp, dp
    implicit none
    private
    public read_mhd_config_from_outfile, init_mhd_data, free_mhd_data, &
           read_mhd_data, save_organized_mhd_data, read_organized_mhd_data, &
           adjust_mhd_data_boundary

    type file_header
        integer, dimension(4) :: nx4
        real(sp) :: time
        real(sp), dimension(4) :: bbox
    end type file_header

    type(file_header) :: fheader

    real(sp), allocatable, dimension(:, :, :, :) :: f_array
    !dir$ attributes align:64 :: f_array

    real(sp), allocatable, dimension(:, :, :, :) :: mhd_data_single_core
    !dir$ attributes align:64 :: mhd_data_single_core

    contains

    !---------------------------------------------------------------------------
    !< Read MHD simulation configuration from the simulation output file
    !< Args:
    !<  filename: the name of one output file
    !<  filename_next_step: the name of the output file of the next time step
    !---------------------------------------------------------------------------
    subroutine read_mhd_config_from_outfile(filename, filename_next_step)
        use mhd_config_module, only: set_mhd_config
        implicit none
        character(*), intent(in) :: filename, filename_next_step
        integer, dimension(4) :: n4
        real(dp) :: t1
        integer :: fh
        real(dp) :: dx, dy, dz, xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: lx, ly, lz, dt_out
        integer :: nx, ny, nz, nxs, nys, nzs
        integer :: topox, topoy, topoz, nvar
        integer :: bcx, bcy, bcz

        fh = 20
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='read')
        read(fh) fheader
        read(fh) n4
        nvar = fheader%nx4(1)
        nx = fheader%nx4(2)
        ny = fheader%nx4(3)
        nxs = n4(3)
        nys = n4(4)
        topox = nx / nxs
        topoy = ny / nys
        xmin = fheader%bbox(1)
        xmax = fheader%bbox(2)
        ymin = fheader%bbox(3)
        ymax = fheader%bbox(4)
        lx = xmax - xmin
        ly = ymax - ymin
        dx = lx / (nx - 1)
        dy = ly / (ny - 1)
        close(fh)

        t1 = fheader%time
        open(unit=fh, file=filename_next_step, access='stream', &
             status='unknown', form='unformatted', action='read')
        read(fh) fheader
        close(fh)
        dt_out = fheader%time - t1

        dz = dx
        zmin = 0.0
        zmax = 1.0
        lz = 1.0
        nz = 1
        nzs = 1
        topoz = 1
        bcx = 0
        bcy = 0
        bcz = 0

        call set_mhd_config(dx, dy, dz, xmin, ymin, zmin, xmax, ymax, zmax, &
            lx, ly, lz, dt_out, nx, ny, nz, nxs, nys, nzs, topox, topoy, topoz, &
            nvar, bcx, bcy, bcz)

    end subroutine read_mhd_config_from_outfile

    !---------------------------------------------------------------------------
    !< Initialize MHD data arrays. Note that we do not read Az, rho and pres.
    !---------------------------------------------------------------------------
    subroutine init_mhd_data
        use mhd_config_module, only: mhd_config
        implicit none
        integer :: nx, ny, nz, nvar, nxs, nys, nzs
        nx = mhd_config%nx
        ny = mhd_config%ny
        nz = mhd_config%nz
        nvar = mhd_config%nvar
        nxs = mhd_config%nxs
        nys = mhd_config%nys
        nzs = mhd_config%nzs

        !< 2D MHD fields are increased by 4 grids along x and y, so 4 ghost
        !< cells are included for calculating the derivatives
        !< vx, vy, vz, pad1, bx, by, bz, btot
        allocate(f_array(8, -1:nx + 2, -1:ny + 2, nz))
        f_array = 0.0

        allocate(mhd_data_single_core(nxs, nys, nzs, nvar))
        mhd_data_single_core = 0.0
    end subroutine init_mhd_data

    !---------------------------------------------------------------------------
    !< Free MHD data arrays
    !---------------------------------------------------------------------------
    subroutine free_mhd_data
        implicit none
        deallocate(f_array)
        deallocate(mhd_data_single_core)
    end subroutine free_mhd_data

    !---------------------------------------------------------------------------
    !< Read MHD data arrays from a file
    !< Args:
    !<  filename: file name to get the data
    !---------------------------------------------------------------------------
    subroutine read_mhd_data(filename)
        use mhd_config_module, only: mhd_config
        implicit none
        character(*), intent(in) :: filename
        integer, dimension(4) :: n4
        integer :: fh, rank, ix, iy, nx1, ny1, ncpus
        fh = 30
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='read')
        read(fh) fheader
        ncpus = mhd_config%topox * mhd_config%topoy
        do rank = 1, ncpus
            read(fh) n4
            ix = n4(1)
            iy = n4(2)
            nx1 = n4(3)
            ny1 = n4(4)
            read(fh) mhd_data_single_core
            if (mhd_config%nvar .eq. 7) then
                f_array(1, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 3, :)
                f_array(2, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 4, :)
                f_array(5, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 5, :)
                f_array(6, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 6, :)
            else
                f_array(1, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 3, :)
                f_array(2, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 4, :)
                f_array(3, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 5, :)
                f_array(5, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 6, :)
                f_array(6, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 7, :)
                f_array(7, ix+1:ix+nx1, iy+1:iy+ny1, :) = mhd_data_single_core(:, :, 8, :)
            endif
        enddo

        !< Calculate the magnitude of the magnetic field
        if (mhd_config%nvar .eq. 7) then
            f_array(8, :, :, :) = sqrt(f_array(5, :, :, :)**2 + f_array(6, :, :, :)**2)
        else
            f_array(8, :, :, :) = sqrt(f_array(5, :, :, :)**2 + &
                                       f_array(6, :, :, :)**2 + &
                                       f_array(7, :, :, :)**2)
        endif
        close(fh)
    end subroutine read_mhd_data

    !---------------------------------------------------------------------------
    !< Adjust MHD data arrays at the boundaries. Currently, only periodic
    !< boundary is carefully implemented for now.
    !---------------------------------------------------------------------------
    subroutine adjust_mhd_data_boundary
        use mhd_config_module, only: mhd_config
        implicit none
        integer :: nx, ny, nz

        nx = mhd_config%nx
        ny = mhd_config%ny
        nz = mhd_config%nz

        if (mhd_config%bcx == 0 .and. mhd_config%bcy == 0) then
            ! 4 borders
            f_array(:, -1:0, :, :) = f_array(:, nx-2:nx-1, :, :)
            f_array(:, nx+1:nx+2, :, :) = f_array(:, 2:3, :, :)
            f_array(:, :, -1:0, :) = f_array(:, :, ny-2:ny-1, :)
            f_array(:, :, ny+1:ny+2, :) = f_array(:, :, 2:3, :)
            ! 4 corners
            f_array(:, -1:0, -1:0, :) = f_array(:, nx-2:nx-1, ny-2:ny-1, :)
            f_array(:, nx+1:nx+2, -1:0, :) = f_array(:, 2:3, ny-2:ny-1, :)
            f_array(:, -1:0, ny+1:ny+2, :) = f_array(:, nx-2:nx-1, 2:3, :)
            f_array(:, nx+1:nx+2, ny+1:ny+2, :) = f_array(:, 2:3, 2:3, :)
        else
            !< Just copy nearby fields for now
            f_array(:, -1, :, :) = f_array(:, 1, :, :)
            f_array(:, 0, :, :) = f_array(:, 1, :, :)
            f_array(:, nx+1, :, :) = f_array(:, nx, :, :)
            f_array(:, nx+2, :, :) = f_array(:, nx, :, :)
            f_array(:, :, -1, :) = f_array(:, :, 1, :)
            f_array(:, :, 0, :) = f_array(:, :, 1, :)
            f_array(:, :, ny+1, :) = f_array(:, :, ny, :)
            f_array(:, :, ny+2, :) = f_array(:, :, ny, :)
        endif
    end subroutine adjust_mhd_data_boundary

    !---------------------------------------------------------------------------
    !< Save re-organized MHD data into a file
    !< Args:
    !<  filename: file name to get the data
    !---------------------------------------------------------------------------
    subroutine save_organized_mhd_data(filename)
        implicit none
        character(*), intent(in) :: filename
        integer :: fh
        fh = 35
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='write')
        write(fh, pos=1) f_array
        close(fh)
    end subroutine save_organized_mhd_data

    !---------------------------------------------------------------------------
    !< Read re-organized MHD data from a file
    !< Args:
    !<  filename: file name to get the data
    !---------------------------------------------------------------------------
    subroutine read_organized_mhd_data(filename)
        implicit none
        character(*), intent(in) :: filename
        integer :: fh
        fh = 35
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='read')
        read(fh, pos=1) f_array
        close(fh)
    end subroutine read_organized_mhd_data

end module mhd_data_sli

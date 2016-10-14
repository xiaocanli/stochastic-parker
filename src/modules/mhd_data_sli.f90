!*******************************************************************************
!< Module of MHD data
!*******************************************************************************
module mhd_data
    use constants, only: fp, dp
    implicit none
    private
    public mhd_config
    public read_mhd_config

    type mhd_configuration
        real(dp) :: dx, dy, xmin, xmax, ymin, ymax, lx, ly
        integer :: nx, ny
    end type mhd_configuration

    type(mhd_configuration) :: mhd_config

    contains

    !---------------------------------------------------------------------------
    !< Read MHD simulation configuration
    !---------------------------------------------------------------------------
    subroutine read_mhd_config
        use read_config, only: get_variable
        implicit none
        real(fp) :: temp
        integer :: fh
        fh = 10
        open(unit=fh, file='config/mhd_config.dat', status='old')
        mhd_config%xmin = get_variable(fh, 'xmin', '=')
        mhd_config%xmax = get_variable(fh, 'xmax', '=')
        mhd_config%ymin = get_variable(fh, 'ymin', '=')
        mhd_config%ymax = get_variable(fh, 'ymax', '=')
        mhd_config%lx = mhd_config%xmax - mhd_config%xmin
        mhd_config%ly = mhd_config%ymax - mhd_config%ymin
        temp = get_variable(fh, 'nx', '=')
        mhd_config%nx = int(temp)
        temp = get_variable(fh, 'ny', '=')
        mhd_config%ny = int(temp)
        mhd_config%dx = mhd_config%lx / mhd_config%nx
        mhd_config%dy = mhd_config%ly / mhd_config%ny
        close(fh)
        ! Echo this information
        print *, "---------------------------------------------------"
        write(*, "(A)") " MHD simulation information."
        write(*, "(A,F7.2,A,F7.2)") " lx, ly = ", &
            mhd_config%lx, ',', mhd_config%ly
        write(*, "(A,I0,A,I0)") " nx, ny = ", &
            mhd_config%nx, ',', mhd_config%ny
        write(*, "(A,F9.6,A,F9.6)") " dx, dy = ", &
            mhd_config%dx, ',', mhd_config%dy
        print *, "---------------------------------------------------"
    end subroutine read_mhd_config
end module mhd_data

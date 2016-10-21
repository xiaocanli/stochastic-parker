!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    implicit none
    private
    public init_particles, free_particles

    type particle_type
        real(dp) :: x, y, p         ! Position and momentum
        real(dp) :: weight, t, dt   ! Particle weight, time and time step
        integer  :: split_times     ! Particle splitting times
        integer  :: count_flag      ! Only count particle when it is 0
    end type particle_type

    type(particle_type), allocatable, dimension(:) :: ptl

    contains

    !---------------------------------------------------------------------------
    !< Initialize particle data
    !---------------------------------------------------------------------------
    subroutine init_particles(nptl_max)
        implicit none
        integer, intent(in) :: nptl_max
        allocate(ptl(nptl_max))
        ptl%x = 0.0
        ptl%y = 0.0
        ptl%p = 0.0
        ptl%weight = 0.0
        ptl%t = 0.0
        ptl%dt = 0.0
        ptl%split_times = 0
        ptl%count_flag = 0
    end subroutine init_particles

    !---------------------------------------------------------------------------
    !< Free particle data
    !---------------------------------------------------------------------------
    subroutine free_particles
        implicit none
        deallocate(ptl)
    end subroutine free_particles

    !---------------------------------------------------------------------------
    !< Move particles using the MHD simulation data as background fields
    !< Args:
    !<  diffusion_type: 1 for constant kappa, 2 for kappa ~ p, 3 for kappa ~ p/B
    !---------------------------------------------------------------------------
    subroutine particle_mover(diffusion_type)
        use mhd_data_sli, only: mhd_config
        implicit none
        integer, intent(in) :: diffusion_type
    end subroutine particle_mover
end module particle_module

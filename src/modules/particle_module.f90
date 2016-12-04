!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    implicit none
    private
    public init_particles, free_particles, inject_particles_spatial_uniform

    type particle_type
        real(dp) :: x, y, p         ! Position and momentum
        real(dp) :: weight, t, dt   ! Particle weight, time and time step
        integer  :: split_times     ! Particle splitting times
        integer  :: count_flag      ! Only count particle when it is 1
    end type particle_type

    type(particle_type), allocatable, dimension(:) :: ptl

    integer :: nptl_current         ! Number of particles currently in the box
    integer :: nptl_max             ! Maximum number of particles allowed

    contains

    !---------------------------------------------------------------------------
    !< Initialize particle data
    !< Args:
    !<  nptl_max_allowed: the maximum number of particles allowd
    !---------------------------------------------------------------------------
    subroutine init_particles(nptl_max_allowed)
        implicit none
        integer, intent(in) :: nptl_max_allowed
        nptl_max = nptl_max_allowed
        allocate(ptl(nptl_max))
        ptl%x = 0.0
        ptl%y = 0.0
        ptl%p = 0.0
        ptl%weight = 0.0
        ptl%t = 0.0
        ptl%dt = 0.0
        ptl%split_times = 0
        ptl%count_flag = 0
        nptl_current = 0     ! No particle initially
    end subroutine init_particles

    !---------------------------------------------------------------------------
    !< Free particle data
    !---------------------------------------------------------------------------
    subroutine free_particles
        implicit none
        deallocate(ptl)
    end subroutine free_particles

    !---------------------------------------------------------------------------
    !< Inject particles which are spatially uniform
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !---------------------------------------------------------------------------
    subroutine inject_particles_spatial_uniform(nptl, dt)
        use mhd_data_sli, only: mhd_config
        use mpi_module, only: mpi_rank, master
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl
        real(dp), intent(in) :: dt
        integer :: i, imod2
        real(dp) :: xmin, ymin, xmax, ymax
        real(dp) :: rands(2)

        xmin = mhd_config%xmin
        xmax = mhd_config%xmax
        ymin = mhd_config%ymin
        ymax = mhd_config%ymax

        do i = 1, nptl
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            ptl(nptl_current)%x = unif_01()*(xmax-xmin) + xmin
            ptl(nptl_current)%y = unif_01()*(ymax-ymin) + ymin
            imod2 = mod(i, 2)
            if (imod2 == 1) rands = two_normals()
            ptl(nptl_current)%p = rands(imod2+1)
            ptl(nptl_current)%weight = 1.0
            ptl(nptl_current)%t = 0.0
            ptl(nptl_current)%dt = dt
            ptl(nptl_current)%split_times = 0
            ptl(nptl_current)%count_flag = 1
        enddo
    end subroutine inject_particles_spatial_uniform

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

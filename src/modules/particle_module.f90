!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    implicit none
    private
    public init_particles, free_particles, inject_particles_spatial_uniform, &
        read_particle_params, particle_mover, remove_particles

    type particle_type
        real(dp) :: x, y, p         !< Position and momentum
        real(dp) :: weight, t, dt   !< Particle weight, time and time step
        integer  :: split_times     !< Particle splitting times
        integer  :: count_flag      !< Only count particle when it is 1
    end type particle_type

    type(particle_type), allocatable, dimension(:) :: ptls
    type(particle_type) :: ptl

    integer :: nptl_current         !< Number of particles currently in the box
    integer :: nptl_max             !< Maximum number of particles allowed
    real(dp) :: leak                !< Leaking particles considering weight

    real(dp) :: kpara0              !< Normalization for kappa parallel
    real(dp) :: kret                !< The ratio of kpara to kperp
    integer :: diffusion_type       !< 1 for constant kappa, 2 for kappa ~ p,
                                    !< 3 for kappa ~ p/B, 4 for kappa ~ p^2/B
    !< These two depends on the simulation normalizations
    real(dp) :: p0  !< the standard deviation of the Gaussian distribution of momentum
    real(dp) :: b0  !< Initial magnetic field strength
    real(dp) :: kpara, kperp, dkxx_dx, dkyy_dy, dkxy_dx, dkxy_dy
    real(dp) :: skperp, skpara_perp

    real(dp) :: dt_min  !< Minimum time step

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
        allocate(ptls(nptl_max))
        ptls%x = 0.0
        ptls%y = 0.0
        ptls%p = 0.0
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = 0
        nptl_current = 0     ! No particle initially
    end subroutine init_particles

    !---------------------------------------------------------------------------
    !< Free particle data
    !---------------------------------------------------------------------------
    subroutine free_particles
        implicit none
        deallocate(ptls)
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
            ptls(nptl_current)%x = unif_01()*(xmax-xmin) + xmin
            ptls(nptl_current)%y = unif_01()*(ymax-ymin) + ymin
            imod2 = mod(i, 2)
            if (imod2 == 1) rands = two_normals()
            ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = 0.0
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = 1
        enddo
    end subroutine inject_particles_spatial_uniform

    !---------------------------------------------------------------------------
    !< Move particles using the MHD simulation data as background fields
    !---------------------------------------------------------------------------
    subroutine particle_mover
        use mhd_data_sli, only: mhd_config, interp_fields, fields, gradf
        implicit none
        real(dp) :: dtf, dxm, dym, xmin, xmax, ymin, ymax
        real(dp) :: t0, px, py, rx, ry, rt, rt1
        real(dp) :: deltax, deltay, deltap
        integer :: i, ix, iy

        dtf = mhd_config%dt_out
        xmin = mhd_config%xmin
        xmax = mhd_config%xmax
        ymin = mhd_config%ymin
        ymax = mhd_config%ymax
        dxm = mhd_config%dx
        dym = mhd_config%dy

        do i = 1, nptl_current
            ptl = ptls(i)
            t0 = ptl%t
            do while ((ptl%t - t0) < dtf .and. ptl%count_flag == 1)
                if(ptl%p < 0.0) exit
                !< Periodic boundary condition
                if (ptl%x > xmax) ptl%x = ptl%x - xmax + xmin
                if (ptl%x < xmin) ptl%x = xmax - (xmin - ptl%x)
                if (ptl%y > ymax) ptl%y = ptl%y - ymax + ymin
                if (ptl%y < ymin) ptl%y = ymax - (ymin - ptl%y)

                px = (ptl%x-xmin) / dxm
                py = (ptl%y-ymin) / dym
                ix = floor(px)
                iy = floor(py)

                rx = px - ix
                ry = py - iy
                rt = (ptl%t - t0) / dtf
                rt1 = 1.0_dp - rt

                call interp_fields(ix, iy, rx, ry, rt)
                call calc_spatial_diffusion_coefficients
                call set_time_step
                call push_particle(rt, deltax, deltay, deltap)
            enddo
           
            if ((ptl%t - t0) > dtf) then
                ptl%x = ptl%x - deltax
                ptl%y = ptl%y - deltay
                ptl%p = ptl%p - deltap
                ptl%t = ptl%t - ptl%dt
                ptl%dt = t0 + dtf - ptl%t

                px = (ptl%x-xmin) / dxm
                py = (ptl%y-ymin) / dym
                ix = floor(px)
                iy = floor(py)

                rx = px - ix
                ry = py - iy
                rt = (ptl%t - t0) / dtf
                rt1 = 1.0_dp - rt

                call interp_fields(ix, iy, rx, ry, rt)
                call calc_spatial_diffusion_coefficients
                call push_particle(rt, deltax, deltay, deltap)
            endif
            !< Periodic boundary condition
            if (ptl%x > xmax) ptl%x = ptl%x - xmax + xmin
            if (ptl%x < xmin) ptl%x = xmax - (xmin - ptl%x)
            if (ptl%y > ymax) ptl%y = ptl%y - ymax + ymin
            if (ptl%y < ymin) ptl%y = ymax - (ymin - ptl%y)

            ptls(i) = ptl
        enddo
    end subroutine particle_mover

    !---------------------------------------------------------------------------
    !< Calculate the spatial diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine calc_spatial_diffusion_coefficients
        use mhd_data_sli, only: fields, gradf
        implicit none
        real(dp) :: pnorm
        real(dp) :: bx, by, b, ib2, ib3, ib4
        real(dp) :: dbx_dx, dby_dx, dbx_dy, dby_dy, db_dx, db_dy

        b = fields%btot
        bx = fields%bx
        by = fields%by
        db_dx = gradf%dbtot_dx
        db_dy = gradf%dbtot_dy
        dbx_dx = gradf%dbx_dx
        dbx_dy = gradf%dbx_dy
        dby_dx = gradf%dby_dx
        dby_dy = gradf%dby_dy
        ib2 = 1.0_dp / b**2
        ib3 = 1.0_dp / b**3
        ib4 = ib2**2
        
        if (diffusion_type == 1) then
            pnorm = 1.0_dp
        else if (diffusion_type == 2) then
            pnorm = ptl%p / p0
        else if (diffusion_type == 3) then
            pnorm = ptl%p * b0 / p0
        else if (diffusion_type == 4) then
            pnorm = ptl%p**2 * b0 / p0**2
        endif

        kpara = kpara0 * pnorm
        kperp = kpara * kret

        if (diffusion_type == 1 .or. diffusion_type == 2) then
            dkxx_dx = 0.0_dp
            dkxy_dy = 0.0_dp
            dkxy_dx = 0.0_dp
            dkyy_dy = 0.0_dp
        else if (diffusion_type == 3 .or. diffusion_type == 4) then
            dkxx_dx = -kperp*db_dx*ib2 + (kperp-kpara)*db_dx*bx**2*ib4 + &
                2.0*(kpara-kperp)*bx*(dbx_dx*b-bx*db_dx)*ib4
            dkyy_dy = -kperp*db_dy*ib2 + (kperp-kpara)*db_dy*by**2*ib4 + &
                2.0*(kpara-kperp)*by*(dby_dy*b-by*db_dy)*ib4
            dkxy_dx = (kperp-kpara)*db_dx*bx*by*ib4 + (kpara-kperp) * &
                ((dbx_dx*by+bx*dby_dx)*ib3 - 2.0*bx*by*db_dx*ib4)
            dkxy_dy = (kperp-kpara)*db_dy*bx*by*ib4 + (kpara-kperp) * &
                ((dbx_dy*by+bx*dby_dy)*ib3 - 2.0*bx*by*db_dy*ib4)
        endif
    end subroutine calc_spatial_diffusion_coefficients

    !---------------------------------------------------------------------------
    !< Read particle parameters including the diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine read_particle_params
        use read_config, only: get_variable
        use mpi_module
        implicit none
        real(fp) :: temp
        integer :: fh
        
        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/conf.dat', status='old')
            p0 = get_variable(fh, 'p0', '=')
            b0 = get_variable(fh, 'b0', '=')
            temp = get_variable(fh, 'diffusion_type', '=')
            diffusion_type = int(temp)
            kpara0 = get_variable(fh, 'kpara0', '=')
            kret = get_variable(fh, 'kret', '=')
            dt_min = get_variable(fh, 'dt_min', '=')
            close(fh)
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A)") " Particle parameters including diffusion coefficients"
            write(*, "(A, E13.6E2)") " The standard deviation of momentum distribution is ", p0
            write(*, "(A, E13.6E2)") " Initial magnetic field strength", b0
            if (diffusion_type == 1) then
                write(*, "(A)") " kpara is a constant"
            else if (diffusion_type == 2) then
                write(*, "(A)") " kpara ~ particle momentum p"
            else if (diffusion_type == 3) then
                write(*, "(A)") " kpara ~ the ratio p / B"
            else if (diffusion_type == 4) then
                write(*, "(A)") " kpara ~ the ratio p^2 / B"
            endif
            write(*, "(A,E13.6E2)") " kpara0 = ", kpara0
            write(*, "(A,E13.6E2)") " kpara / kperp = ", kret
            write(*, "(A,E13.6E2)") " Minimum time step = ", dt_min
            print *, "---------------------------------------------------"
        endif
        call MPI_BCAST(diffusion_type, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p0, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(b0, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kpara0, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kret, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_min, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
    end subroutine read_particle_params

    !---------------------------------------------------------------------------
    !< Determine the time step.
    !---------------------------------------------------------------------------
    subroutine set_time_step
        use mhd_data_sli, only: fields, mhd_config
        implicit none
        real(dp) :: tmp30, tmp40, bx, by, b, vx, vy, dxm, dym, dt1, dt2
        skperp = dsqrt(2.0*kperp)
        skpara_perp = dsqrt(2.0*(kpara-kperp))

        vx = fields%vx
        vy = fields%vy
        bx = fields%bx
        by = fields%by
        b = fields%btot
        dxm = mhd_config%dx
        dym = mhd_config%dy

        tmp30 = skperp + skpara_perp * abs(bx/b)
        tmp40 = abs(vx + dkxx_dx + dkxy_dy)
        if (tmp40 .ne. 0.0d0) then
            dt1 = min(dxm/(80.0*tmp40), (tmp30/(tmp40))**2) * 0.5
        else
            dt1 = dt_min
        endif
        tmp30 = skperp + skpara_perp * abs(by/b)
        tmp40 = abs(vy + dkxy_dx + dkyy_dy)
        if (tmp40 .ne. 0.0d0) then
            dt2 = min(dym/(80.0*tmp40), (tmp30/(tmp40))**2) * 0.5
        else
            dt2 = dt_min
        endif
        ptl%dt = min(dt1, dt2)
        if (ptl%dt .lt. dt_min) then
            ptl%dt = dt_min
        endif
    end subroutine set_time_step

    !---------------------------------------------------------------------------
    !< Push particle for a single step
    !< Args:
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  deltax, deltay, deltap: the change of x, y and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle(rt, deltax, deltay, deltap)
        use mhd_data_sli, only: fields, gradf, mhd_config, interp_fields
        use random_number_generator, only: unif_01
        implicit none
        real(dp), intent(in) :: rt
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp) :: xtmp, ytmp, ptmp
        real(dp) :: sdt, dvx_dx, dvy_dy
        real(dp) :: bx, by, b, vx, vy, px, py, rx, ry, rt1
        real(dp) :: bx1, by1, b1
        real(dp) :: xmin, ymin, dxm, dym, skperp1, skpara_perp1
        real(dp) :: ran1, ran2, ran3, sqrt3
        integer :: ix, iy

        b = fields%btot
        bx = fields%bx
        by = fields%by
        vx = fields%vx
        vy = fields%vy
        dvx_dx = gradf%dvx_dx
        dvy_dy = gradf%dvy_dy
        xmin = mhd_config%xmin
        ymin = mhd_config%ymin
        dxm = mhd_config%dx
        dym = mhd_config%dy

        sdt = dsqrt(ptl%dt)
        xtmp = ptl%x + (vx+dkxx_dx+dkxy_dy)*ptl%dt + (skperp+skpara_perp*bx/b)*sdt
        ytmp = ptl%y + (vy+dkxy_dx+dkyy_dy)*ptl%dt + (skperp+skpara_perp*by/b)*sdt
        ptmp = ptl%p - ptl%p*(dvx_dx+dvy_dy)*ptl%dt/3.0d0

        if (ptmp .lt. 0.0_dp) then
            ptl%count_flag = 0
            leak = leak + ptl%weight
            return
        endif

        px = (xtmp - xmin) / dxm
        py = (ytmp - ymin) / dym
        ix = floor(px)
        iy = floor(py)
        rx = px - ix
        ry = py - iy
        rt1 = 1.0_dp - rt
        call interp_fields(ix, iy, rx, ry, rt)
        call calc_spatial_diffusion_coefficients

        !< Magnetic field at the predicted position
        b1 = fields%btot
        bx1 = fields%bx
        by1 = fields%by

        skperp1 = dsqrt(2.0*kperp)
        skpara_perp1 = dsqrt(2.0*(kpara-kperp))

        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3

        deltax = (vx+dkxx_dx+dkxy_dy)*ptl%dt + ran1*skperp*sdt + &
                 ran3*skpara_perp*sdt*bx/b + &
                 (skperp1-skperp)*(ran1**2-1.0)*sdt/2.0 + &
                 (skpara_perp1*bx1/b1-skpara_perp*bx/b)*(ran3**2-1.0)*sdt/2.0
        deltay = (vy+dkxy_dx+dkyy_dy)*ptl%dt + ran2*skperp*sdt + &
                 ran3*skpara_perp*sdt*by/b + &
                 (skperp1-skperp)*(ran2**2-1.0)*sdt/2.0 + &
                 (skpara_perp1*by1/b1-skpara_perp*by/b)*(ran3**2-1.0)*sdt/2.0
        deltap = -ptl%p * (dvx_dx+dvy_dy) * ptl%dt / 3.0d0
        ptl%x = ptl%x + deltax 
        ptl%y = ptl%y + deltay 
        ptl%p = ptl%p + deltap
        ptl%t = ptl%t + ptl%dt
    end subroutine push_particle

    !---------------------------------------------------------------------------
    !< Remove particles from simulation if their count_flags are 0.
    !---------------------------------------------------------------------------
    subroutine remove_particles
        implicit none
        integer :: i
        i = 1
        do while (i <= nptl_current)
            if (ptls(i)%count_flag == 0) then
                !< Switch current with the last particle
                ptl = ptls(nptl_current)
                ptls(nptl_current) = ptls(i)
                ptls(i) = ptl
                nptl_current = nptl_current - 1
            else
                i = i + 1
            endif
        enddo
    end subroutine remove_particles
end module particle_module

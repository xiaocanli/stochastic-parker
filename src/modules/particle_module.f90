!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    implicit none
    private
    public init_particles, free_particles, inject_particles_spatial_uniform, &
        read_particle_params, particle_mover, remove_particles, split_particle, &
        init_particle_distributions, clean_particle_distributions, &
        free_particle_distributions, distributions_diagnostics, quick_check

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
    integer :: nptl_new             !< Number of particles from splitting
    real(dp) :: leak                !< Leaking particles considering weight

    real(dp) :: kpara0              !< Normalization for kappa parallel
    real(dp) :: kret                !< The ratio of kpara to kperp
    integer :: diffusion_type       !< 1 for constant kappa, 2 for kappa ~ p,
                                    !< 3 for kappa ~ p/B, 4 for kappa ~ p^2/B
    !< These two depends on the simulation normalizations
    real(dp) :: p0  !< the standard deviation of the Gaussian distribution of momentum
    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum
    real(dp) :: b0    !< Initial magnetic field strength
    real(dp) :: kpara, kperp, dkxx_dx, dkyy_dy, dkxy_dx, dkxy_dy
    real(dp) :: skperp, skpara_perp

    real(dp) :: dt_min  !< Minimum time step

    !< Dimensions of particle distributions
    integer :: nx, ny, npp, nreduce
    real(dp) :: dx, dy, pmin_log, pmax_log, dp_log

    real(dp), allocatable, dimension(:, :) :: f0, f1, f2, f3
    real(dp), allocatable, dimension(:, :) :: f0_sum, f1_sum, f2_sum, f3_sum
    real(dp), allocatable, dimension(:) :: fp0, fp0_sum
    real(dp), allocatable, dimension(:, :) :: fp1, fp1_sum
    real(dp), allocatable, dimension(:, :) :: fp2, fp2_sum
    real(dp), allocatable, dimension(:) :: fx0, fx0_sum
    real(dp), allocatable, dimension(:) :: fy0, fy0_sum
    real(dp), allocatable, dimension(:) :: parray

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
    !<  dist_flag: momentum distribution flag. 0 for Maxwellian, 1 for delta.
    !---------------------------------------------------------------------------
    subroutine inject_particles_spatial_uniform(nptl, dt, dist_flag)
        use mhd_data_sli, only: mhd_config
        use mpi_module, only: mpi_rank, master
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag
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
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals()
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = 0.0
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = 1
        enddo
        leak = 0.0_dp
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
                ix = floor(px) + 1
                iy = floor(py) + 1

                rx = px + 1 - ix
                ry = py + 1 - iy
                rt = (ptl%t - t0) / dtf
                rt1 = 1.0_dp - rt

                call interp_fields(ix, iy, rx, ry, rt)
                call calc_spatial_diffusion_coefficients
                call set_time_step(t0, dtf)
                call push_particle(rt, deltax, deltay, deltap)
            enddo

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
        use mhd_data_sli, only: mhd_config
        use mpi_module
        implicit none
        real(fp) :: temp
        integer :: fh
        
        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/conf.dat', status='old')
            b0 = get_variable(fh, 'b0', '=')
            p0 = get_variable(fh, 'p0', '=')
            pmin = get_variable(fh, 'pmin', '=')
            pmax = get_variable(fh, 'pmax', '=')
            temp = get_variable(fh, 'diffusion_type', '=')
            diffusion_type = int(temp)
            kpara0 = get_variable(fh, 'kpara0', '=')
            kret = get_variable(fh, 'kret', '=')
            dt_min = get_variable(fh, 'dt_min', '=')
            temp = get_variable(fh, 'nreduce', '=')
            nreduce = int(temp)
            nx = mhd_config%nx / nreduce
            ny = mhd_config%ny / nreduce
            temp = get_variable(fh, 'npp', '=')
            npp = int(temp)
            close(fh)
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A)") " Particle parameters including diffusion coefficients"
            write(*, "(A, E13.6E2)") " The standard deviation of momentum distribution is ", p0
            write(*, "(A, E13.6E2, E13.6E2)") &
                " Minimum and maximum particle momentum", pmin, pmax
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
            write(*, "(A,I0,A,I0)") " Dimensions of spatial distributions = ", &
                nx, " ", ny
            write(*, "(A,I0)") " Dimensions of momentum distributions = ", npp
            print *, "---------------------------------------------------"
        endif
        call MPI_BCAST(diffusion_type, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p0, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pmin, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pmax, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(b0, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kpara0, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kret, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_min, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nreduce, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nx, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ny, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(npp, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    end subroutine read_particle_params

    !---------------------------------------------------------------------------
    !< Determine the time step.
    !< Args;
    !<  t0: the initial time for current particle
    !<  dtf: the time interval of the MHD fields
    !---------------------------------------------------------------------------
    subroutine set_time_step(t0, dtf)
        use mhd_data_sli, only: fields, mhd_config
        implicit none
        real(dp), intent(in) :: t0, dtf
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
            dt1 = min(dxm/(80.0*tmp40), tmp30/(tmp40**2)) * 0.5
        else
            dt1 = dt_min
        endif
        tmp30 = skperp + skpara_perp * abs(by/b)
        tmp40 = abs(vy + dkxy_dx + dkyy_dy)
        if (tmp40 .ne. 0.0d0) then
            dt2 = min(dym/(80.0*tmp40), tmp30/(tmp40**2)) * 0.5
        else
            dt2 = dt_min
        endif
        ptl%dt = min(dt1, dt2)

        !< Make sure the time step is not too small
        if (ptl%dt .lt. dt_min) then
            ptl%dt = dt_min
        endif

        !< Make sure the time step is not too large
        if ((ptl%t + ptl%dt - t0) > dtf) then
            ptl%dt = t0 + dtf - ptl%t
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
        real(dp) :: xmin, ymin, xmax, ymax, dxm, dym, skperp1, skpara_perp1
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: lx, ly
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
        xmax = mhd_config%xmax
        ymax = mhd_config%ymax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        lx = mhd_config%lx
        ly = mhd_config%ly

        sdt = dsqrt(ptl%dt)
        xtmp = ptl%x + (vx+dkxx_dx+dkxy_dy)*ptl%dt + (skperp+skpara_perp*bx/b)*sdt
        ytmp = ptl%y + (vy+dkxy_dx+dkyy_dy)*ptl%dt + (skperp+skpara_perp*by/b)*sdt
        ptmp = ptl%p - ptl%p*(dvx_dx+dvy_dy)*ptl%dt/3.0d0

        if (ptmp .lt. 0.0_dp) then
            ptl%count_flag = 0
            leak = leak + ptl%weight
            return
        endif

        !< Periodic boundary condition
        if (xtmp > xmax) xtmp = xtmp - xmax + xmin
        if (xtmp < xmin) xtmp = xmax - (xmin - xtmp)
        if (ytmp > ymax) ytmp = ytmp - ymax + ymin
        if (ytmp < ymin) ytmp = ymax - (ymin - ytmp)

        px = (xtmp - xmin) / dxm
        py = (ytmp - ymin) / dym
        ix = floor(px) + 1
        iy = floor(py) + 1
        rx = px + 1 - ix
        ry = py + 1 - iy
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

    !---------------------------------------------------------------------------
    !< When a particle get to a certain energy, split it to two particles 
    !< to increase statistics at high energy
    !---------------------------------------------------------------------------
    subroutine split_particle
        implicit none
        integer :: i, nptl
        real(dp) :: p_threshold
        type(particle_type) :: ptl_new
        nptl = nptl_current
        do i = 1, nptl
            ptl = ptls(i)
            p_threshold = (1+2.72**ptl%split_times)*p0
            if (ptl%p > p_threshold .and. ptl%p <= pmax) then
                nptl_current = nptl_current + 1
                if (nptl_current > nptl_max) then
                    nptl_current = nptl_max
                    return
                endif
                nptl_new = nptl_new + 1
                ptl%weight = 0.5**(1.0 + ptl%split_times)
                ptl%split_times = ptl%split_times + 1
                ptls(nptl_current) = ptl
                ptls(i) = ptl
            endif
        enddo
    end subroutine split_particle

    !---------------------------------------------------------------------------
    !< Quick diagnostics of the particle number information
    !< Args:
    !<  iframe: the time frame
    !<  if_create_file: whether to create a file
    !---------------------------------------------------------------------------
    subroutine quick_check(iframe, if_create_file)
        use mpi_module, only: mpi_rank, master
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: if_create_file
        real(dp) :: pdt_min, pdt_max, pdt_avg, ntot
        integer :: i
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif
        if (mpi_rank == master) then
            ntot = 0.0_dp
            pdt_min = 1.0_dp
            pdt_max = 0.0_dp
            pdt_avg = 0.0_dp
            do i = 1, nptl_current
                ntot = ntot + ptls(i)%weight
                if (ptls(i)%dt < pdt_min) pdt_min = ptls(i)%dt
                if (ptls(i)%dt > pdt_max) pdt_max = ptls(i)%dt
                pdt_avg = pdt_avg + ptls(i)%dt
            enddo
            pdt_avg = pdt_avg / nptl_current
            if (if_create_file) then
                open (17, file='data/quick.dat', status='unknown')
                write(17, "(A,A)") "iframe, nptl_current, nptl_new, ntot, ", &
                    "leak, pdt_min, pdt_max, pdt_avg"
            else
                open (17, file='data/quick.dat', status="old", &
                    position="append", action="write")
            endif
            write(17, "(I4.4,A,I6.6,A,I6.6,A,5E13.6E2)") &
                iframe, ' ', nptl_current, ' ', nptl_new, ' ', ntot, &
                leak, pdt_min, pdt_max, pdt_avg
            close(17)
        endif 
    end subroutine quick_check

    !---------------------------------------------------------------------------
    !< Initialize the particle distributions
    !---------------------------------------------------------------------------
    subroutine init_particle_distributions
        use mhd_data_sli, only: mhd_config
        use mpi_module, only: mpi_rank, master
        implicit none
        integer :: i
        allocate(f0(nx, ny))
        allocate(f1(nx, ny))
        allocate(f2(nx, ny))
        allocate(f3(nx, ny))
        allocate(fp0(npp))
        allocate(fp1(npp, nx))
        allocate(fp2(npp, ny))
        allocate(fx0(nx))
        allocate(fy0(ny))
        allocate(f0_sum(nx, ny))
        allocate(f1_sum(nx, ny))
        allocate(f2_sum(nx, ny))
        allocate(f3_sum(nx, ny))
        allocate(fp0_sum(npp))
        allocate(fp1_sum(npp, nx))
        allocate(fp2_sum(npp, ny))
        allocate(fx0_sum(nx))
        allocate(fy0_sum(ny))
        call clean_particle_distributions

        !< Intervals for distributions
        dx = mhd_config%lx / (mhd_config%nx / dble(nreduce) - 1)
        dy = mhd_config%ly / (mhd_config%ny / dble(nreduce) - 1)
        pmin_log = log10(pmin)
        pmax_log = log10(pmax)
        dp_log = (pmax_log - pmin_log) / (npp - 1)

        if (mpi_rank == master) then
            allocate(parray(npp))
            parray = 0.0_dp
            do i = 1, npp
                parray(i) = 10**(pmin_log + (i-1) * dp_log)
            enddo
        endif
    end subroutine init_particle_distributions

    !---------------------------------------------------------------------------
    !< Set particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_particle_distributions
        implicit none
        f0 = 0.0_dp
        f1 = 0.0_dp
        f2 = 0.0_dp
        f3 = 0.0_dp
        fp0 = 0.0_dp
        fp1 = 0.0_dp
        fp2 = 0.0_dp
        fx0 = 0.0_dp
        fy0 = 0.0_dp
        f0_sum = 0.0_dp
        f1_sum = 0.0_dp
        f2_sum = 0.0_dp
        f3_sum = 0.0_dp
        fp0_sum = 0.0_dp
        fp1_sum = 0.0_dp
        fp2_sum = 0.0_dp
        fx0_sum = 0.0_dp
        fy0_sum = 0.0_dp
    end subroutine clean_particle_distributions

    !---------------------------------------------------------------------------
    !< Free particle distributions
    !---------------------------------------------------------------------------
    subroutine free_particle_distributions
        use mpi_module, only: mpi_rank, master
        implicit none
        deallocate(f0, f1, f2, f3, fp0, fp1, fp2, fx0, fy0)
        deallocate(f0_sum, f1_sum, f2_sum, f3_sum)
        deallocate(fp0_sum, fp1_sum, fp2_sum, fx0_sum, fy0_sum)
        if (mpi_rank == master) then
            deallocate(parray)
        endif
    end subroutine free_particle_distributions

    !---------------------------------------------------------------------------
    !< Accumulate particle distributions
    !---------------------------------------------------------------------------
    subroutine calc_particle_distributions
        use mpi_module
        use mhd_data_sli, only: mhd_config
        implicit none
        integer :: i, ix, iy
        real(dp) :: weight, p, ip, xmin, xmax, ymin, ymax

        xmin = mhd_config%xmin
        xmax = mhd_config%xmax
        ymin = mhd_config%ymin
        ymax = mhd_config%ymax
        
        call clean_particle_distributions

        do i = 1, nptl_current
            ptl = ptls(i)
            ix = ceiling((ptl%x - xmin)/dx)
            iy = ceiling((ptl%y - ymin)/dy)
            if (ix < 1) ix = 1
            if (ix > nx) ix = nx
            if (iy < 1) iy = 1
            if (iy > ny) iy = ny
            p = ptl%p
            weight = ptl%weight
            fx0(ix) = fx0(ix) + weight
            fy0(iy) = fy0(iy) + weight

            !< Different momentum band
            if (p > 0.5*p0 .and. p <= 1.5*p0) then
                f0(ix,iy) = f0(ix,iy) + weight
            endif
            if (p > 1.5*p0 .and. p <= 2.5*p0) then
                f1(ix,iy) = f1(ix,iy) + weight
            endif
            if (p > 2.5*p0 .and. p <= 3.5*p0) then
                f2(ix,iy) = f2(ix,iy) + weight
            endif
            if (p > 3.5*p0 .and. p <= 4.5*p0) then
                f3(ix,iy) = f3(ix,iy) + weight
            endif

            if (p > pmin .and. p <= pmax) then
                ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
                fp0(ip) = fp0(ip) + weight
                fp1(ip, ix) = fp1(ip, ix) + weight
                fp2(ip, iy) = fp2(ip, iy) + weight
            endif
        enddo
        
        call MPI_REDUCE(fx0, fx0_sum, nx, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fy0, fy0_sum, ny, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f0, f0_sum, nx*ny, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f1, f1_sum, nx*ny, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f2, f2_sum, nx*ny, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f3, f3_sum, nx*ny, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp0, fp0_sum, npp, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp1, fp1_sum, npp*nx, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp2, fp2_sum, npp*ny, MPI_DOUBLE, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
    end subroutine calc_particle_distributions

    !---------------------------------------------------------------------------
    !< Diagnostics of the particle distributions
    !< Args:
    !<  iframe: time frame index
    !---------------------------------------------------------------------------
    subroutine distributions_diagnostics(iframe)
        use mpi_module
        implicit none
        integer, intent(in) :: iframe
        integer :: ix, iy, ip, fh
        character(len=4) :: ctime
        character(len=64) :: fname
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif

        call calc_particle_distributions

        write (ctime,'(i4.4)') iframe
        if (mpi_rank .eq. 0) then
            fh = 13
            fname = 'data/fx-'//ctime//'.dat'
            open(fh, file=trim(fname), status='unknown')
            do ix = 1, nx
                write (fh, "(E13.6E2)") fx0_sum(ix)
            enddo
            close(fh)

            fh = 14
            fname = 'data/fy-'//ctime//'.dat'
            open(fh, file=trim(fname), status='unknown')
            do iy = 1, ny
                write (fh, "(E13.6E2)") fy0_sum(iy)
            enddo
            close(fh)

            fh = 15
            fname = 'data/fp-'//ctime//'.dat'
            open(fh, file=trim(fname), status='unknown')
            do ip = 1, npp
                write (fh, "(2E13.6E2)") parray(ip), fp0_sum(ip)
            enddo
            close(fh)

            fh = 16
            fname = 'data/fpx-'//ctime//'.dat'
            open(fh, file=trim(fname), status='unknown')
            do ix = 1, nx
                do ip = 1, npp
                    write (fh, "(2E13.6E2)") fp1_sum(ip, ix)
                enddo
            enddo
            close(fh)

            fh = 17
            fname = 'data/fpy-'//ctime//'.dat'
            open(fh, file=trim(fname), status='unknown')
            do iy = 1, ny
                do ip = 1, npp
                    write (fh, "(2E13.6E2)") fp2_sum(ip, iy)
                enddo
            enddo
            close(fh)

            fh = 18
            fname = 'data/fxy-'//ctime//'.dat'
            open(fh, file=trim(fname), status='unknown')
            do iy = 1, ny
                do ix = 1, nx
                    write (fh, "(4E13.6E2)") f0_sum(ix, iy), f1_sum(ix, iy),&
                        f2_sum(ix, iy), f3_sum(ix, iy)
                enddo
            enddo
            close(fh)
        endif
    end subroutine distributions_diagnostics
end module particle_module

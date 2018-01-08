!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    implicit none
    private
    save
    public init_particles, free_particles, inject_particles_spatial_uniform, &
        read_particle_params, particle_mover, remove_particles, split_particle, &
        init_particle_distributions, clean_particle_distributions, &
        free_particle_distributions, distributions_diagnostics, quick_check, &
        set_particle_datatype_mpi, free_particle_datatype_mpi, &
        select_particles_tracking, init_particle_tracking, &
        free_particle_tracking, init_tracked_particle_points, &
        free_tracked_particle_points, negative_particle_tags, &
        save_tracked_particle_points, inject_particles_at_shock

    type particle_type
        real(dp) :: x, y, p         !< Position and momentum
        real(dp) :: weight, t, dt   !< Particle weight, time and time step
        integer  :: split_times     !< Particle splitting times
        integer  :: count_flag      !< Only count particle when it is 1
        integer  :: tag             !< Particle tag
        integer  :: nsteps_tracking !< Total particle tracking steps
    end type particle_type

    integer :: particle_datatype_mpi
    type(particle_type), allocatable, dimension(:) :: ptls
    type(particle_type), allocatable, dimension(:, :) :: senders
    type(particle_type), allocatable, dimension(:, :) :: recvers
    integer, dimension(4) :: nsenders, nrecvers
    !dir$ attributes align:64 :: ptls
    type(particle_type) :: ptl, ptl1
    
    ! Particle tracking
    integer, allocatable, dimension(:) :: nsteps_tracked_ptls, noffsets_tracked_ptls
    integer, allocatable, dimension(:) :: tags_selected_ptls
    type(particle_type), allocatable, dimension(:) :: ptl_traj_points
    integer :: nsteps_tracked_tot   !< Total tracking steps for all tracked particles

    integer :: nptl_current         !< Number of particles currently in the box
    integer :: nptl_old             !< Number of particles without receivers
    integer :: nptl_max             !< Maximum number of particles allowed
    integer :: nptl_new             !< Number of particles from splitting
    real(dp) :: leak                !< Leaking particles considering weight

    real(dp) :: kpara0              !< Normalization for kappa parallel
    real(dp) :: kret                !< The ratio of kpara to kperp
    integer :: momentum_dependency  !< kappa dependency on particle momentum
    integer :: mag_dependency       !< kappa dependency on magnetic field
    real(dp) :: pindex              !< power index for the momentum dependency 
    real(dp) :: p0    !< the standard deviation of the Gaussian distribution of momentum
    real(dp) :: b0    !< Initial magnetic field strength
    real(dp) :: kpara, kperp, dkxx_dx, dkyy_dy, dkxy_dx, dkxy_dy
    real(dp) :: skperp, skpara_perp

    real(dp) :: dt_min  !< Minimum time step

    !< Parameters for particle distributions
    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum
    real(dp) :: dx, dy, pmin_log, pmax_log, dp_log
    integer :: nx, ny, npp, nreduce

    real(dp), allocatable, dimension(:, :) :: f0, f1, f2, f3
    real(dp), allocatable, dimension(:, :) :: f0_sum, f1_sum, f2_sum, f3_sum
    real(dp), allocatable, dimension(:) :: fp0, fp0_sum
    real(dp), allocatable, dimension(:, :) :: fp1, fp1_sum
    real(dp), allocatable, dimension(:, :) :: fp2, fp2_sum
    real(dp), allocatable, dimension(:) :: parray

    contains

    !---------------------------------------------------------------------------
    !< Initialize particle data
    !< Args:
    !<  nptl_max_allowed: the maximum number of particles allowed
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
        ptls%tag = 0
        ptls%nsteps_tracking = 0
        nptl_current = 0     ! No particle initially

        !< Particles crossing domain boundaries
        allocate(senders(nptl_max / 100, 4))
        allocate(recvers(nptl_max / 100, 4))
        senders%x = 0.0
        senders%y = 0.0
        senders%p = 0.0
        senders%weight = 0.0
        senders%t = 0.0
        senders%dt = 0.0
        senders%split_times = 0
        senders%count_flag = 0
        senders%tag = 0
        senders%nsteps_tracking = 0

        recvers%x = 0.0
        recvers%y = 0.0
        recvers%p = 0.0
        recvers%weight = 0.0
        recvers%t = 0.0
        recvers%dt = 0.0
        recvers%split_times = 0
        recvers%count_flag = 0
        recvers%tag = 0
        recvers%nsteps_tracking = 0

        nsenders = 0
        nrecvers = 0
    end subroutine init_particles

    !---------------------------------------------------------------------------
    !< Free particle data
    !---------------------------------------------------------------------------
    subroutine free_particles
        implicit none
        deallocate(ptls, senders, recvers)
    end subroutine free_particles

    !---------------------------------------------------------------------------
    !< Set MPI datatype for particle type
    !---------------------------------------------------------------------------
    subroutine set_particle_datatype_mpi
        use mpi_module
        implicit none
        integer :: oldtypes(0:1), blockcounts(0:1)
        integer :: offsets(0:1), extent
        ! Setup description of the 8 MPI_DOUBLE fields.
        offsets(0) = 0
        oldtypes(0) = MPI_DOUBLE_PRECISION
        blockcounts(0) = 6
        ! Setup description of the 7 MPI_INTEGER fields.
        call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
        offsets(1) = blockcounts(0) * extent
        oldtypes(1) = MPI_INTEGER
        blockcounts(1) = 4
        ! Define structured type and commit it. 
        call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
            particle_datatype_mpi, ierr)
        call MPI_TYPE_COMMIT(particle_datatype_mpi, ierr)
    end subroutine set_particle_datatype_mpi

    !---------------------------------------------------------------------------
    !< Free MPI datatype for particle type
    !---------------------------------------------------------------------------
    subroutine free_particle_datatype_mpi
        use mpi_module
        implicit none
        call MPI_TYPE_FREE(particle_datatype_mpi, ierr)
    end subroutine free_particle_datatype_mpi

    !---------------------------------------------------------------------------
    !< Inject particles which are spatially uniform
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: momentum distribution flag. 0 for Maxwellian, 1 for delta.
    !---------------------------------------------------------------------------
    subroutine inject_particles_spatial_uniform(nptl, dt, dist_flag)
        use simulation_setup_module, only: fconfig
        use mpi_module, only: mpi_rank, master
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag
        real(dp), intent(in) :: dt
        integer :: i, imod2
        real(dp) :: xmin, ymin, xmax, ymax
        real(dp) :: rands(2)

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax

        nptl_current = 0
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
            ptls(nptl_current)%tag = nptl_current
        enddo
        leak = 0.0_dp

        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles"
        endif
    end subroutine inject_particles_spatial_uniform

    !---------------------------------------------------------------------------
    !< Inject particles at shock
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: momentum distribution flag. 0 for Maxwellian, 1 for delta.
    !<  ct_mhd: MHD simulation time frame
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_shock(nptl, dt, dist_flag, ct_mhd)
        use simulation_setup_module, only: fconfig
        use mpi_module, only: mpi_rank, master
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_shock_location
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt
        integer :: i, imod2, iy, ix
        real(dp) :: xmin, ymin, xmax, ymax, dpy, shock_xpos
        real(dp), dimension(2) :: rands

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax

        do i = 1, nptl
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            ptls(nptl_current)%y = unif_01()*(ymax-ymin) + ymin
            dpy = ptls(nptl_current)%y / mhd_config%dy
            iy = floor(dpy)
            !< We assume that particles are inject at the earlier shock location
            shock_xpos = interp_shock_location(iy, dpy - iy, 0.0d0) + 2  ! Two ghost cells
            ptls(nptl_current)%x = shock_xpos * (xmax - xmin) / fconfig%nxg
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals()
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = 1
            ptls(nptl_current)%tag = nptl_current
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles at the shock"
        endif
    end subroutine inject_particles_at_shock

    !---------------------------------------------------------------------------
    !< Particle mover in one cycle
    !< Args:
    !<  t0: the starting time
    !<  track_particle_flag: whether to track particles, 0 for no, 1 for yes
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !---------------------------------------------------------------------------
    subroutine particle_mover_one_cycle(t0, track_particle_flag, nptl_selected, &
            nsteps_interval)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields
        use mpi_module
        implicit none
        real(dp), intent(in) :: t0
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        real(dp) :: dtf, dxm, dym, xmin, xmax, ymin, ymax
        real(dp) :: px, py, rx, ry, rt, rt1
        real(dp) :: deltax, deltay, deltap
        integer :: i, ix, iy, j, tracking_step, offset

        dtf = mhd_config%dt_out
        dxm = mhd_config%dx
        dym = mhd_config%dy
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax

        j = 0

        do i = nptl_old + 1, nptl_current
            ptl = ptls(i)
            do while ((ptl%t - t0) < dtf .and. ptl%count_flag /= 0)
                ! if (mpi_rank == 12) then
                !     print*, ptl%t, t0
                ! endif
                if (ptl%p < 0.0) then
                    ptl%count_flag = 0
                    leak = leak + ptl%weight
                    exit
                else
                    call particle_boundary_condition(ptl%x, ptl%y, xmin, xmax, ymin, ymax)
                    if (ptl%count_flag == 0) then
                        exit
                    endif
                endif
                call particle_boundary_condition(ptl%x, ptl%y, xmin, xmax, ymin, ymax)
                if (ptl%count_flag == 0) then
                    exit
                endif

                j = j + 1
                ! if (nptl_current == nptl_old + 1) then
                !     print*, mpi_rank, ptl%t, ptl%dt, ptl%tag
                ! endif

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
                
                ! Number of particle tracking steps
                ptl%nsteps_tracking = ptl%nsteps_tracking + 1

                ! Track particles
                if (track_particle_flag == 1) then
                    if (ptl%tag < 0 .and. &
                        mod(ptl%nsteps_tracking, nsteps_interval) == 0) then
                        tracking_step = ptl%nsteps_tracking / nsteps_interval
                        offset = noffsets_tracked_ptls(-ptl%tag)
                        ptl_traj_points(offset + tracking_step + 1) = ptl
                    endif
                endif
            enddo

            ! if (nptl_current == nptl_old + 1) then
            !     print*, mpi_rank, j
            ! endif

            if ((ptl%t - t0) > dtf .and. ptl%count_flag /= 0) then
                ptl%x = ptl%x - deltax
                ptl%y = ptl%y - deltay
                ptl%p = ptl%p - deltap
                ptl%t = ptl%t - ptl%dt
                ptl%dt = t0 + dtf - ptl%t
                ptl%nsteps_tracking = ptl%nsteps_tracking - 1

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
                call push_particle(rt, deltax, deltay, deltap)
                ptl%nsteps_tracking = ptl%nsteps_tracking + 1

                ! Track particles
                if (track_particle_flag == 1) then
                    if (ptl%tag < 0 .and. &
                        mod(ptl%nsteps_tracking, nsteps_interval) == 0) then
                        tracking_step = ptl%nsteps_tracking / nsteps_interval
                        offset = noffsets_tracked_ptls(-ptl%tag)
                        ptl_traj_points(offset + tracking_step + 1) = ptl
                    endif
                endif

                if (ptl%p < 0.0) then
                    ptl%count_flag = 0
                    leak = leak + ptl%weight
                else
                    call particle_boundary_condition(ptl%x, ptl%y, xmin, xmax, ymin, ymax)
                endif
            endif

            ptls(i) = ptl

        enddo
    end subroutine particle_mover_one_cycle

    !---------------------------------------------------------------------------
    !< Move particles using the MHD simulation data as background fields
    !< Args:
    !<  track_particle_flag: whether to track particles, 0 for no, 1 for yes
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !---------------------------------------------------------------------------
    subroutine particle_mover(track_particle_flag, nptl_selected, nsteps_interval)
        use simulation_setup_module, only: fconfig
        use mpi_module
        implicit none
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        integer :: i, local_flag, global_flag, ncycle
        logical :: all_particles_in_box
        real(dp) :: t0
        all_particles_in_box = .false.
        nptl_old = 0
        nptl_new = 0

        t0 = ptls(1)%t

        ncycle = 0
        local_flag = 0
        global_flag = 0
        
        do while (.not. all_particles_in_box)
            ncycle = ncycle + 1
            nsenders = 0
            nrecvers = 0
            if (nptl_old < nptl_current) then
                call particle_mover_one_cycle(t0, track_particle_flag, &
                    nptl_selected, nsteps_interval)
            endif
            call remove_particles
            call send_recv_particles
            call add_neighbor_particles
            if (sum(nrecvers) > 0) then
                local_flag = sum(nrecvers)
            else
                local_flag = 0
            endif
            call MPI_ALLREDUCE(local_flag, global_flag, 1, MPI_INTEGER, &
                MPI_SUM, MPI_COMM_WORLD, ierr)
            if (global_flag > 0) then
                all_particles_in_box = .false.
                ! if (global_flag == 1 .and. (mpi_rank == 1 .or. mpi_rank == 2)) then
                !     print*, mpi_rank, nsenders(1:2), nrecvers(1:2)
                !     print*, mpi_rank, recvers(1, 1:2)%x
                !     print*, mpi_rank, recvers(1, 1:2)%tag
                !     print*, mpi_rank, recvers(1, 1:2)%dt
                !     print*, mpi_rank, recvers(1, 1:2)%t
                !     print*, mpi_rank, senders(1, 1:2)%x
                !     print*, mpi_rank, senders(1, 1:2)%tag
                !     print*, mpi_rank, senders(1, 1:2)%dt
                !     print*, mpi_rank, senders(1, 1:2)%t
                !     print*, mpi_rank, ptl%t, ptl%dt
                !     print*, '--------'
                ! endif
                ! if (mpi_rank == master) then
                !     print*, ncycle, global_flag
                ! endif
            else
                all_particles_in_box = .true.
            endif
        enddo
        if (mpi_rank == master) then
            write(*, "(A, I0)") "Number of cycles: ", ncycle
        endif
    end subroutine particle_mover

    !---------------------------------------------------------------------------
    !< Particle boundary conditions
    !< Args:
    !<  x, y: particle positions
    !<  xmin, xmax: min and max along the x-direction
    !<  ymin, ymax: min and max along the y-direction
    !---------------------------------------------------------------------------
    subroutine particle_boundary_condition(x, y, xmin, xmax, ymin, ymax)
        use simulation_setup_module, only: neighbors, mpi_ix, mpi_iy, &
            mpi_sizex, mpi_sizey
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        real(dp), intent(inout) :: x, y

        if (x < xmin) then
            if (neighbors(1) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0
            else if (neighbors(1) == mpi_rank) then
                x = x - xmin + xmax
            else if (neighbors(1) == mpi_rank + mpi_sizex - 1) then
                x = x - mhd_config%xmin + mhd_config%xmax
                nsenders(1) = nsenders(1) + 1
                senders(nsenders(1), 1) = ptl
                ptl%count_flag = 0
            else
                nsenders(1) = nsenders(1) + 1
                senders(nsenders(1), 1) = ptl
                ptl%count_flag = 0
            endif
        else if (x > xmax) then
            if (neighbors(2) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0 !< remove particle
            else if (neighbors(2) == mpi_rank) then
                x = x - xmax + xmin
            else if (neighbors(2) == mpi_iy * mpi_sizex) then !< simulation boundary
                x = x - mhd_config%xmax + mhd_config%xmin
                nsenders(2) = nsenders(2) + 1
                senders(nsenders(2), 2) = ptl
                ptl%count_flag = 0
            else
                nsenders(2) = nsenders(2) + 1
                senders(nsenders(2), 2) = ptl
                ptl%count_flag = 0
            endif
        else if (y < ymin) then
            if (neighbors(3) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0
            else if (neighbors(3) == mpi_rank) then
                y = y - ymin + ymax
            else if (neighbors(3) == mpi_rank + (mpi_sizey - 1) * mpi_sizex) then
                y = y - mhd_config%ymin + mhd_config%ymax
                nsenders(3) = nsenders(3) + 1
                senders(nsenders(3), 3) = ptl
                ptl%count_flag = 0
            else
                nsenders(3) = nsenders(3) + 1
                senders(nsenders(3), 3) = ptl
                ptl%count_flag = 0
            endif
        else if (y > ymax) then
            if (neighbors(4) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0
            else if (neighbors(4) == mpi_rank) then
                y = y - ymax + ymin
            else if (neighbors(4) == mpi_ix) then
                y = y - mhd_config%ymax + mhd_config%ymin
                nsenders(4) = nsenders(4) + 1
                senders(nsenders(4), 4) = ptl
                ptl%count_flag = 0
            else
                nsenders(4) = nsenders(4) + 1
                senders(nsenders(4), 4) = ptl
                ptl%count_flag = 0
            endif
        endif
    end subroutine particle_boundary_condition

    !---------------------------------------------------------------------------
    !< Send and receiver particles from one neighbor
    !< Args:
    !<  send_id, recv_id: sender and receiver ID
    !<  mpi_direc: MPI rank along one direction
    !---------------------------------------------------------------------------
    subroutine send_recv_particle_one_neighbor(send_id, recv_id, mpi_direc)
        use mpi_module
        use simulation_setup_module, only: neighbors
        implicit none
        integer, intent(in) :: send_id, recv_id, mpi_direc
        integer :: nsend, nrecv
        nrecv = 0
        nsend = nsenders(send_id)
        if (neighbors(send_id) /= mpi_rank) then
            call MPI_SEND(nsend, 1, MPI_INTEGER, neighbors(send_id), &
                mpi_rank, MPI_COMM_WORLD, ierr)
        endif
        if (neighbors(recv_id) /= mpi_rank) then
            call MPI_RECV(nrecv, 1, MPI_INTEGER, &
                neighbors(recv_id), neighbors(recv_id), MPI_COMM_WORLD, status, ierr)
        endif
        nrecvers(recv_id) = nrecv
        !< This assumes MPI size along this direction is even
        if (mpi_direc / 2 == 0) then
            if (nsend > 0) then
                call MPI_SEND(senders(1:nsend, send_id), nsend, &
                    particle_datatype_mpi, neighbors(send_id), mpi_rank, &
                    MPI_COMM_WORLD, ierr)
            endif
            if (nrecv > 0) then
                call MPI_RECV(recvers(1:nrecv, recv_id), nrecv, &
                    particle_datatype_mpi, neighbors(recv_id), &
                    neighbors(recv_id), MPI_COMM_WORLD, status, ierr)
            endif
        else
            if (nrecv > 0) then
                call MPI_RECV(recvers(1:nrecv, recv_id), nrecv, &
                    particle_datatype_mpi, neighbors(recv_id), &
                    neighbors(recv_id), MPI_COMM_WORLD, status, ierr)
            endif
            if (nsend > 0) then
                call MPI_SEND(senders(1:nsend, send_id), nsend, &
                    particle_datatype_mpi, neighbors(send_id), mpi_rank, &
                    MPI_COMM_WORLD, ierr)
            endif
        endif
    end subroutine send_recv_particle_one_neighbor

    !---------------------------------------------------------------------------
    !< Send particles to neighbors
    !---------------------------------------------------------------------------
    subroutine send_recv_particles
        use mpi_module
        use simulation_setup_module, only: neighbors, mpi_ix, mpi_iy
        implicit none
        call send_recv_particle_one_neighbor(1, 2, mpi_ix) !< Right -> Left
        call send_recv_particle_one_neighbor(2, 1, mpi_ix) !< Left -> Right
        call send_recv_particle_one_neighbor(3, 4, mpi_iy) !< Top -> Bottom
        call send_recv_particle_one_neighbor(4, 3, mpi_iy) !< Bottom -> Top
    end subroutine send_recv_particles

    !---------------------------------------------------------------------------
    !< Calculate the spatial diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine calc_spatial_diffusion_coefficients
        use mhd_data_parallel, only: fields, gradf
        implicit none
        real(dp) :: pnorm
        real(dp) :: bx, by, b, ib1, ib2, ib3, ib4
        real(dp) :: dbx_dx, dby_dx, dbx_dy, dby_dy, db_dx, db_dy

        bx = fields(5)
        by = fields(6)
        b = fields(8)
        dbx_dx = gradf(4)
        dbx_dy = gradf(5)
        dby_dx = gradf(7)
        dby_dy = gradf(8)
        db_dx = gradf(13)
        db_dy = gradf(14)
        ib1 = 1.0_dp / b
        ib2 = ib1 * ib1
        ib3 = ib1 * ib2
        ib4 = ib2 * ib2
        
        pnorm = 1.0_dp
        if (mag_dependency == 1) then
            pnorm = pnorm * b0
        endif
        if (momentum_dependency == 1) then
            pnorm = pnorm * (ptl%p / p0)**pindex
        endif

        kpara = kpara0 * pnorm
        kperp = kpara * kret

        if (mag_dependency == 1) then
            dkxx_dx = -kperp*db_dx*ib2 + (kperp-kpara)*db_dx*bx*bx*ib4 + &
                2.0*(kpara-kperp)*bx*(dbx_dx*b-bx*db_dx)*ib4
            dkyy_dy = -kperp*db_dy*ib2 + (kperp-kpara)*db_dy*by*by*ib4 + &
                2.0*(kpara-kperp)*by*(dby_dy*b-by*db_dy)*ib4
            dkxy_dx = (kperp-kpara)*db_dx*bx*by*ib4 + (kpara-kperp) * &
                ((dbx_dx*by+bx*dby_dx)*ib3 - 2.0*bx*by*db_dx*ib4)
            dkxy_dy = (kperp-kpara)*db_dy*bx*by*ib4 + (kpara-kperp) * &
                ((dbx_dy*by+bx*dby_dy)*ib3 - 2.0*bx*by*db_dy*ib4)
        else
            dkxx_dx = 2.0*(kpara-kperp)*bx*(dbx_dx*b-bx*db_dx)*ib3
            dkyy_dy = 2.0*(kpara-kperp)*by*(dby_dy*b-by*db_dy)*ib3
            dkxy_dx = (kpara-kperp) * ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
            dkxy_dy = (kpara-kperp) * ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
        endif
    end subroutine calc_spatial_diffusion_coefficients

    !---------------------------------------------------------------------------
    !< Read particle parameters including the diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine read_particle_params
        use read_config, only: get_variable
        use simulation_setup_module, only: fconfig
        use simulation_setup_module, only: mpi_sizex, mpi_sizey
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
            temp = get_variable(fh, 'momentum_dependency', '=')
            momentum_dependency = int(temp)
            if (momentum_dependency == 1) then
                pindex = get_variable(fh, 'pindex', '=')
            endif
            temp = get_variable(fh, 'mag_dependency', '=')
            mag_dependency = int(temp)
            kpara0 = get_variable(fh, 'kpara0', '=')
            kret = get_variable(fh, 'kret', '=')
            dt_min = get_variable(fh, 'dt_min', '=')
            temp = get_variable(fh, 'nreduce', '=')
            nreduce = int(temp)
            nx = fconfig%nx / nreduce
            ny = fconfig%ny / nreduce
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
            write(*, "(A,I0)") " kappa dependency on p: ", momentum_dependency
            if (momentum_dependency == 1) then
                write(*, "(A,E13.6E2)") " kappa ~ p^", pindex
            else
                write(*, "(A)") " kappa doesn't depend on particle momentum"
            endif
            if (mag_dependency == 1) then
                write(*, "(A)") " kappa ~ 1/B"
            else
                write(*, "(A)") " kappa doesn't depend on B"
            endif
            write(*, "(A,E13.6E2)") " kpara0 = ", kpara0
            write(*, "(A,E13.6E2)") " kpara / kperp = ", kret
            write(*, "(A,E13.6E2)") " Minimum time step = ", dt_min
            write(*, "(A,I0,A,I0)") " Dimensions of spatial distributions = ", &
                nx, " ", ny
            write(*, "(A,I0)") " Dimensions of momentum distributions = ", npp
            print *, "---------------------------------------------------"
        endif
        call MPI_BCAST(momentum_dependency, 1, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mag_dependency, 1, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pindex, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p0, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pmin, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pmax, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(b0, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kpara0, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kret, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_min, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nreduce, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nx, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ny, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(npp, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        if (nx * nreduce /= fconfig%nx .or. ny * nreduce /= fconfig%ny) then
            if (mpi_rank == master) then 
                write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
            endif
            call MPI_FINALIZE(ierr)
            stop
        endif
    end subroutine read_particle_params

    !---------------------------------------------------------------------------
    !< Determine the time step.
    !< Args;
    !<  t0: the initial time for current particle
    !<  dtf: the time interval of the MHD fields
    !---------------------------------------------------------------------------
    subroutine set_time_step(t0, dtf)
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: fields
        implicit none
        real(dp), intent(in) :: t0, dtf
        real(dp) :: tmp30, tmp40, bx, by, b, vx, vy, dxm, dym, dt1, dt2
        skperp = dsqrt(2.0*kperp)
        skpara_perp = dsqrt(2.0*(kpara-kperp))

        vx = fields(1)
        vy = fields(2)
        bx = fields(5)
        by = fields(6)
        b = fields(8)
        dxm = mhd_config%dx
        dym = mhd_config%dy

        tmp30 = skperp + skpara_perp * abs(bx/b)
        tmp40 = abs(vx + dkxx_dx + dkxy_dy)
        if (tmp40 .ne. 0.0d0) then
            ! dt1 = min(dxm/(80.0*tmp40), (tmp30/tmp40)**2)
            dt1 = min(dxm/tmp40, (dxm/tmp30) * (dxm/tmp30))
        else
            dt1 = dt_min
        endif
        tmp30 = skperp + skpara_perp * abs(by/b)
        tmp40 = abs(vy + dkxy_dx + dkyy_dy)
        if (tmp40 .ne. 0.0d0) then
            ! dt2 = min(dym/(80.0*tmp40), (tmp30/tmp40)**2)
            dt2 = min(dym/tmp40, (dym/tmp30) * (dym/tmp30))
        else
            dt2 = dt_min
        endif
        ptl%dt = min(dt1, dt2)

        !< Make sure the time step is not too large. Adding dt_min to make
        !< sure to exit the while where this routine is called
        if ((ptl%t + ptl%dt - t0) > dtf) then
            ptl%dt = t0 + dtf - ptl%t + dt_min
        endif

        !< Make sure the time step is not too small
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
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: fields, gradf, interp_fields
        use random_number_generator, only: unif_01, two_normals
        implicit none
        real(dp), intent(in) :: rt
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp) :: xtmp, ytmp
        real(dp) :: sdt, dvx_dx, dvy_dy
        real(dp) :: bx, by, b, vx, vy, px, py, rx, ry, rt1
        real(dp) :: bx1, by1, b1
        real(dp) :: xmin, ymin, xmax, ymax, dxm, dym, skperp1, skpara_perp1
        reaL(dp) :: xmin1, ymin1, xmax1, ymax1, dxmh, dymh
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rands(2)
        integer :: ix, iy

        vx = fields(1)
        vy = fields(2)
        bx = fields(5)
        by = fields(6)
        b = fields(8)
        dvx_dx = gradf(1)
        dvy_dy = gradf(2)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dxmh = 0.5 * dxm
        dymh = 0.5 * dym

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        ymin1 = ymin - dymh
        xmax1 = xmax + dxmh
        ymax1 = ymax + dymh

        deltax = 0.0
        deltay = 0.0
        deltap = 0.0

        sdt = dsqrt(ptl%dt)
        xtmp = ptl%x + (vx+dkxx_dx+dkxy_dy)*ptl%dt + (skperp+skpara_perp*bx/b)*sdt
        ytmp = ptl%y + (vy+dkxy_dx+dkyy_dy)*ptl%dt + (skperp+skpara_perp*by/b)*sdt

        !< Make sure the point is still in the box
        do while (xtmp < xmin1 .or. xtmp > xmax1 .or. ytmp < ymin1 .or. ytmp > ymax1)
            ptl%dt = ptl%dt * 0.5
            sdt = dsqrt(ptl%dt)
            xtmp = ptl%x + (vx+dkxx_dx+dkxy_dy)*ptl%dt + (skperp+skpara_perp*bx/b)*sdt
            ytmp = ptl%y + (vy+dkxy_dx+dkyy_dy)*ptl%dt + (skperp+skpara_perp*by/b)*sdt
        enddo

        px = (xtmp - xmin) / dxm
        py = (ytmp - ymin) / dym
        ix = floor(px) + 1
        iy = floor(py) + 1
        rx = px + 1 - ix
        ry = py + 1 - iy
        rt1 = 1.0_dp - rt
        call interp_fields(ix, iy, rx, ry, rt)

        ! Before dkxx_dx, dkxy_dy, dkxy_dx, dkyy_dy are updated
        deltax = (vx+dkxx_dx+dkxy_dy)*ptl%dt
        deltay = (vy+dkxy_dx+dkyy_dy)*ptl%dt

        call calc_spatial_diffusion_coefficients

        !< Magnetic field at the predicted position
        bx1 = fields(5)
        by1 = fields(6)
        b1 = fields(8)

        skperp1 = dsqrt(2.0*kperp)
        skpara_perp1 = dsqrt(2.0*(kpara-kperp))

        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3

        ! rands = two_normals()
        ! ran1 = rands(1)
        ! ran2 = rands(2)
        ! rands = two_normals()
        ! ran3 = rands(1)

        deltax = deltax + ran1*skperp*sdt + ran3*skpara_perp*sdt*bx/b + &
                 (skperp1-skperp)*(ran1*ran1-1.0)*sdt/2.0 + &
                 (skpara_perp1*bx1/b1-skpara_perp*bx/b)*(ran3*ran3-1.0)*sdt/2.0
        deltay = deltay + ran2*skperp*sdt + ran3*skpara_perp*sdt*by/b + &
                 (skperp1-skperp)*(ran2*ran2-1.0)*sdt/2.0 + &
                 (skpara_perp1*by1/b1-skpara_perp*by/b)*(ran3*ran3-1.0)*sdt/2.0
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
                ptl1 = ptls(nptl_current)
                ptls(nptl_current) = ptls(i)
                ptls(i) = ptl1
                nptl_current = nptl_current - 1
            else
                i = i + 1
            endif
        enddo
    end subroutine remove_particles

    !---------------------------------------------------------------------------
    !< Add particles from neighbors
    !---------------------------------------------------------------------------
    subroutine add_neighbor_particles
        implicit none
        integer :: i, j, nrecv
        nptl_old = nptl_current
        do i = 1, 4
            nrecv = nrecvers(i)
            if (nrecv > 0) then
                ptls(nptl_current+1:nptl_current+nrecv) = recvers(1:nrecv, i)
                nptl_current = nptl_current + nrecv
            endif
        enddo
    end subroutine add_neighbor_particles

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
                ptls(nptl_current)%tag = nptl_current
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
        use mhd_config_module, only: mhd_config
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
        allocate(f0_sum(nx, ny))
        allocate(f1_sum(nx, ny))
        allocate(f2_sum(nx, ny))
        allocate(f3_sum(nx, ny))
        allocate(fp0_sum(npp))
        allocate(fp1_sum(npp, nx))
        allocate(fp2_sum(npp, ny))
        call clean_particle_distributions

        !< Intervals for distributions
        dx = mhd_config%lx / (mhd_config%nx / dble(nreduce) - 1)
        dy = mhd_config%ly / (mhd_config%ny / dble(nreduce) - 1)
        pmin_log = log10(pmin)
        pmax_log = log10(pmax)
        dp_log = (pmax_log - pmin_log) / (npp - 1)

        ! if (mpi_rank == master) then
            allocate(parray(npp))
            parray = 0.0_dp
            do i = 1, npp
                parray(i) = 10**(pmin_log + (i-1) * dp_log)
            enddo
        ! endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished Initializing particle distributions."
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
        f0_sum = 0.0_dp
        f1_sum = 0.0_dp
        f2_sum = 0.0_dp
        f3_sum = 0.0_dp
        fp0_sum = 0.0_dp
        fp1_sum = 0.0_dp
        fp2_sum = 0.0_dp
    end subroutine clean_particle_distributions

    !---------------------------------------------------------------------------
    !< Free particle distributions
    !---------------------------------------------------------------------------
    subroutine free_particle_distributions
        use mpi_module, only: mpi_rank, master
        implicit none
        deallocate(f0, f1, f2, f3, fp0, fp1, fp2)
        deallocate(f0_sum, f1_sum, f2_sum, f3_sum)
        deallocate(fp0_sum, fp1_sum, fp2_sum)
        ! if (mpi_rank == master) then
            deallocate(parray)
        ! endif
    end subroutine free_particle_distributions

    !---------------------------------------------------------------------------
    !< Accumulate particle distributions
    !---------------------------------------------------------------------------
    subroutine calc_particle_distributions
        use mpi_module
        use simulation_setup_module, only: fconfig
        implicit none
        integer :: i, ix, iy, ip
        real(dp) :: weight, p, xmin, xmax, ymin, ymax

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        
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
        
        call MPI_REDUCE(f0, f0_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f1, f1_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f2, f2_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(f3, f3_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp0, fp0_sum, npp, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp1, fp1_sum, npp*nx, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp2, fp2_sum, npp*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
    end subroutine calc_particle_distributions

    !---------------------------------------------------------------------------
    !< Diagnostics of the particle distributions
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine distributions_diagnostics(iframe, file_path)
        use constants, only: fp, dp
        use mpi_module
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer :: ix, iy, ip, fh, pos1
        character(len=4) :: ctime, mrank
        character(len=128) :: fname
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif

        call calc_particle_distributions

        write (ctime,'(i4.4)') iframe
        write (mrank,'(i4.4)') mpi_rank
        if (mpi_rank .eq. 0) then
            fh = 15
            ! fname = trim(file_path)//'fp-'//ctime//'_'//mrank//'.dat'
            fname = trim(file_path)//'fp-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fp0_sum
            close(fh)

            ! fh = 16
            ! fname = 'data/fpx-'//ctime//'_'//mrank//'.dat'
            ! open(fh, file=trim(fname), access='stream', status='unknown', &
            !      form='unformatted', action='write')
            ! write(fh, pos=1) fp1
            ! close(fh)

            ! fh = 17
            ! fname = 'data/fpy-'//ctime//'_'//mrank//'.dat'
            ! open(fh, file=trim(fname), access='stream', status='unknown', &
            !      form='unformatted', action='write')
            ! write(fh, pos=1) fp2
            ! close(fh)

            fh = 18
            ! fname = trim(file_path)//'fxy-'//ctime//'_'//mrank//'.dat'
            fname = trim(file_path)//'fxy-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) f0_sum
            pos1 = nx * ny * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) f1_sum
            pos1 = pos1 + nx * ny * sizeof(1.0_dp)
            write(fh, pos=pos1) f2_sum
            pos1 = pos1 + nx * ny * sizeof(1.0_dp)
            write(fh, pos=pos1) f3_sum
            close(fh)
        endif
    end subroutine distributions_diagnostics

    !---------------------------------------------------------------------------
    !< quicksort particles w.r.t their energy/momentum
    !---------------------------------------------------------------------------
    recursive subroutine quicksort_particle(ptls, first, last)
        implicit none
        type(particle_type), dimension(:) :: ptls
        type(particle_type) :: ptl_tmp
        real(dp) :: p_pivot
        integer :: first, last
        integer :: i, j

        p_pivot = ptls((first+last) / 2)%p
        i = first
        j = last
        do
            do while (ptls(i)%p < p_pivot)
                i = i + 1
            end do
            do while (p_pivot < ptls(j)%p)
                j = j - 1
            end do
            if (i >= j) exit
            ptl_tmp = ptls(i)
            ptls(i) = ptls(j)
            ptls(j) = ptl_tmp
            i = i + 1
            j = j - 1
        end do
        if (first < i-1) call quicksort_particle(ptls, first, i-1)
        if (j+1 < last)  call quicksort_particle(ptls, j+1, last)
    end subroutine quicksort_particle

    !---------------------------------------------------------------------------
    !< Initialize particle tracking
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine init_particle_tracking(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        allocate(nsteps_tracked_ptls(nptl_selected))
        allocate(noffsets_tracked_ptls(nptl_selected))
        allocate(tags_selected_ptls(nptl_selected))
        nsteps_tracked_ptls = 0
        noffsets_tracked_ptls = 0
        tags_selected_ptls = 0
    end subroutine init_particle_tracking

    !---------------------------------------------------------------------------
    !< Free particle tracking
    !---------------------------------------------------------------------------
    subroutine free_particle_tracking
        implicit none
        deallocate(nsteps_tracked_ptls, noffsets_tracked_ptls)
        deallocate(tags_selected_ptls)
    end subroutine free_particle_tracking

    !---------------------------------------------------------------------------
    !< Select particles for particle tracking
    !< Args:
    !<  nptl: number of initial particles
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !---------------------------------------------------------------------------
    subroutine select_particles_tracking(nptl, nptl_selected, nsteps_interval)
        implicit none
        integer, intent(in) :: nptl, nptl_selected, nsteps_interval
        integer :: iptl, i
        call quicksort_particle(ptls, 1, nptl_current)
        
        !< Select high-energy particles
        iptl = nptl_current-nptl_selected+1
        tags_selected_ptls = ptls(iptl:nptl_current)%tag
        nsteps_tracked_ptls = ptls(iptl:nptl_current)%nsteps_tracking

        !< Only part of the trajectory points are saved
        nsteps_tracked_ptls = (nsteps_tracked_ptls + nsteps_interval - 1) / nsteps_interval
        nsteps_tracked_ptls = nsteps_tracked_ptls + 1  ! Include the initial point

        do i = 2, nptl_selected
            noffsets_tracked_ptls(i) = noffsets_tracked_ptls(i-1) + &
                nsteps_tracked_ptls(i-1)
        enddo
        nsteps_tracked_tot = noffsets_tracked_ptls(nptl_selected) + &
                             nsteps_tracked_ptls(nptl_selected)

        !< Set particles to their initial values, but switch the tags for
        !< for high-energy particles
        ptls%x = 0.0
        ptls%y = 0.0
        ptls%p = 0.0
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = 0
        ptls%tag = 0
        ptls%nsteps_tracking = 0
    end subroutine select_particles_tracking

    !---------------------------------------------------------------------------
    !< Make the flags of tracked particles negative, so they can be easily
    !< identified.
    !< that can be easily tracked.
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine negative_particle_tags(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        integer :: i
        do i = 1, nptl_selected
            ptls(tags_selected_ptls(i))%tag = -i
        enddo
    end subroutine negative_particle_tags

    !---------------------------------------------------------------------------
    !< Initialize tracked particle points
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine init_tracked_particle_points(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        integer :: i, offset, tag
        allocate(ptl_traj_points(nsteps_tracked_tot))
        ptl_traj_points%x = 0.0
        ptl_traj_points%y = 0.0
        ptl_traj_points%p = 0.0
        ptl_traj_points%weight = 0.0
        ptl_traj_points%t = 0.0
        ptl_traj_points%dt = 0.0
        ptl_traj_points%split_times = 0
        ptl_traj_points%count_flag = 0
        ptl_traj_points%tag = 0
        ptl_traj_points%nsteps_tracking = 0

        !< Save the initial information of tracked particles
        do i = 1, nptl_selected
            offset = noffsets_tracked_ptls(i)
            tag = tags_selected_ptls(i)
            ptl_traj_points(offset + 1) = ptls(tag)
        enddo
    end subroutine init_tracked_particle_points

    !---------------------------------------------------------------------------
    !< Free tracked particle points
    !---------------------------------------------------------------------------
    subroutine free_tracked_particle_points
        implicit none
        deallocate(ptl_traj_points)
    end subroutine free_tracked_particle_points

    !---------------------------------------------------------------------------
    !< Save tracked particle points
    !< Args:
    !<  nptl_selected: number of selected particles
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_tracked_particle_points(nptl_selected, file_path)
        use mpi_module
        implicit none
        integer, intent(in) :: nptl_selected
        character(*), intent(in) :: file_path
        character(len=4) :: mrank
        character(len=64) :: fname
        integer :: fh, offset
        write (mrank,'(i4.4)') mpi_rank
        fh = 41
        fname = trim(file_path)//'tracked_particle_points_'//mrank//'.dat'
        open(unit=fh, file=fname, access='stream', status='unknown', &
            form='unformatted', action='write')     
        write(fh, pos=1) nptl_selected
        write(fh, pos=5) nsteps_tracked_ptls
        offset = (nptl_selected + 1) * sizeof(fp)
        write(fh, pos=offset+1) ptl_traj_points
        close(fh)
    end subroutine save_tracked_particle_points
end module particle_module

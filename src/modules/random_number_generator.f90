!*******************************************************************************
!< Module for generating random numbers for applications using MPI
!< The random number generator is Mersenne Twister. We use multiple stream
!< Mersenne Twister PRNG to generate random numbers for each MPI process.
!< http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html
!*******************************************************************************
module random_number_generator
    use constants, only: dp
    use mt_stream
    implicit none
    private
    public init_prng, delete_prng, unif_01, two_normals
    type(mt_state) :: mts
    type(mt_state), allocatable, dimension(:) :: mtss
    integer :: iseeda(4) = (/ int(Z'123'), int(Z'234'), int(Z'345'), int(Z'456') /)

    contains

    !---------------------------------------------------------------------------
    !< Initialize pseudo random number generator
    !< Args:
    !<  mpi_rank: MPI rank
    !<  nthreads: number of threads on current rank
    !---------------------------------------------------------------------------
    subroutine init_prng(mpi_rank, nthreads)
        implicit none
        integer, intent(in) :: mpi_rank, nthreads
        integer :: i, offset
        call set_mt19937
        call new(mts)
    !    call init(mts,iseed)  ! init by scalar
        call init(mts,iseeda)  ! init by array
        ! if (mpi_rank == master) call print(mts)

        allocate(mtss(nthreads))
        do i = 1, nthreads
            offset = mpi_rank * nthreads + i - 1
            call create_stream(mts, mtss(i), offset)
        enddo
        ! if (mpi_rank == master) call print(mts1)
    end subroutine init_prng

    !---------------------------------------------------------------------------
    !< Delete pseudo random number generator
    !---------------------------------------------------------------------------
    subroutine delete_prng
        implicit none
        integer :: i, nthreads
        call delete(mts)
        nthreads = ubound(mtss, 1) 
        do i = 1, nthreads
            call delete(mtss(i))
        enddo
        deallocate(mtss)
    end subroutine delete_prng

    !---------------------------------------------------------------------------
    !< Uniform number in [0, 1] with double precision
    !< Args:
    !<  thread: the thread starting from 0
    !---------------------------------------------------------------------------
    real(dp) function unif_01(thread)
        implicit none
        integer, intent(in) :: thread
        unif_01 = genrand_double1(mtss(thread+1))
    end function unif_01

    !---------------------------------------------------------------------------
    !< Generate two normal random variates with mean 0 and variance 1.
    !< http://en.wikipedia.org/wiki/Marsaglia_polar_method
    !< Args:
    !<  thread: the thread starting from 0
    !---------------------------------------------------------------------------
    function two_normals(thread) result(rands)
        implicit none
        integer, intent(in) :: thread
        real(dp) :: rands(2), sum_sq
        do
            rands(1) = 2 * unif_01(thread) - 1
            rands(2) = 2 * unif_01(thread) - 1
            sum_sq = sum(rands**2)
            if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
        end do
        rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
    end function two_normals
end module random_number_generator

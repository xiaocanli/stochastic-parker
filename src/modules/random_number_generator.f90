!*******************************************************************************
!< Module for generating random numbers for applications using MPI
!< The random number generator is Mersenne Twister. We use multiple stream
!< Mersenne Twister PRNG to generate random numbers for each MPI process.
!< http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html
!*******************************************************************************
module random_number_generator
    use constants, only: dp
    use mpi_module
    use mt_stream
    implicit none
    private
    public mts, mts1, init_prng, delete_prng, unif_01, two_normals
    type(mt_state) :: mts, mts1
    integer :: iseed = 123456789
    integer :: iseeda(4) = (/ Z'123', Z'234', Z'345', Z'456' /)

    contains

    !---------------------------------------------------------------------------
    !< Initialize pseudo random number generator
    !---------------------------------------------------------------------------
    subroutine init_prng
        implicit none
        call set_mt19937
        call new(mts)
    !    call init(mts,iseed)  ! init by scalar
        call init(mts,iseeda)  ! init by array
        if (mpi_rank == master) call print(mts)

        call create_stream(mts, mts1, mpi_rank)
        if (mpi_rank == master) call print(mts1)
    end subroutine init_prng

    !---------------------------------------------------------------------------
    !< Delete pseudo random number generator
    !---------------------------------------------------------------------------
    subroutine delete_prng
        implicit none

        call delete(mts)
        call delete(mts1)
    end subroutine delete_prng

    !---------------------------------------------------------------------------
    !< Uniform number in [0, 1] with double precision
    !---------------------------------------------------------------------------
    real(dp) function unif_01()
        implicit none
        unif_01 = genrand_double1(mts1)
    end function unif_01

    !---------------------------------------------------------------------------
    !< Generate two normal random variates with mean 0 and variance 1.
    !< http://en.wikipedia.org/wiki/Marsaglia_polar_method
    !---------------------------------------------------------------------------
    function two_normals() result(rands)
        implicit none
        real(dp) :: rands(2), sum_sq
        do
            rands(1) = 2 * unif_01() - 1
            rands(2) = 2 * unif_01() - 1
            sum_sq = sum(rands**2)
            if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
        end do
        rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
    end function two_normals
end module random_number_generator

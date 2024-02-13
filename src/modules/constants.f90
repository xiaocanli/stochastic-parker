!*******************************************************************************
!< Module of the constants used in this code.
!*******************************************************************************
module constants
    use iso_fortran_env
    implicit none
    private
    public sp, dp, qp, i1, i2, i4, i8, pi
    integer, parameter :: sp        = real32
    integer, parameter :: dp        = real64
    integer, parameter :: qp        = real128
    integer, parameter :: i1        = int8
    integer, parameter :: i2        = int16
    integer, parameter :: i4        = int32
    integer, parameter :: i8        = int64

    real(dp), parameter :: pi=4*atan(1.0_dp)
end module constants



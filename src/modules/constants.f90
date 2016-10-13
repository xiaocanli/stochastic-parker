!*******************************************************************************
!< Module of the constants used in this code.
!*******************************************************************************
module constants
    implicit none
    private
    public dp, fp, delta, qi, qe
    integer, parameter :: dp=kind(0.d0)
    integer, parameter :: fp=kind(0.0)
    real(dp), parameter :: delta = 1.0E-15      ! Tiny number
    real(fp), parameter :: qi = 1.0, qe = -1.0  ! Charges
end module constants



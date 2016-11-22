!*******************************************************************************
!< Module for MPI info
!*******************************************************************************
module mpi_module
    implicit none
    include "mpif.h"
    integer :: mpi_rank, mpi_size, ierr
    integer, parameter :: master = 0

    integer :: ierror, ierror2, err_length
    integer :: status(MPI_STATUS_SIZE)
    character(len=256) :: err_msg
end module mpi_module



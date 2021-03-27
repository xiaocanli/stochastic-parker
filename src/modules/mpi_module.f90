!*******************************************************************************
!< Module for MPI info
!< MPI_COMM_WORLD is divided into n_sub_comm. The sub_mpi_rank in each
!< mpi_sub_comm is grouped into mpi_cross_comm.
!*******************************************************************************
module mpi_module
    implicit none
    include "mpif.h"
    integer :: mpi_rank, mpi_size, ierr
    integer :: mpi_sub_comm, mpi_cross_comm
    integer :: mpi_sub_rank, mpi_sub_size
    integer :: mpi_cross_rank, mpi_cross_size
    integer, parameter :: master = 0

    integer :: ierror, ierror2, err_length
    integer :: status(MPI_STATUS_SIZE)
    character(len=256) :: err_msg
end module mpi_module

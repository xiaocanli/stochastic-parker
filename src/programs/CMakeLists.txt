#include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# # For a large number of source files you can create it in a simpler way
# # using file() function:
# file(GLOB analysis_SOURCES *.f90)

# Create an executable files.
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    add_executable(stochastic-mhd.exec stochastic-mhd.f90
        $ENV{MT_STREAM_GCC}/mt_stream.o $ENV{MT_STREAM_GCC}/f_jump_ahead_coeff/gf2xe.o
        $ENV{MT_STREAM_GCC}/f_jump_ahead_coeff/f_get_coeff.o)
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    add_executable(stochastic-mhd.exec stochastic-mhd.f90
        $ENV{MT_STREAM}/mt_stream.o $ENV{MT_STREAM}/f_jump_ahead_coeff/gf2xe.o
        $ENV{MT_STREAM}/f_jump_ahead_coeff/f_get_coeff.o)
elseif (Fortran_COMPILER_NAME MATCHES "ftn*")
    add_executable(stochastic-mhd.exec stochastic-mhd.f90
        $ENV{MT_STREAM_HSW}/mt_stream.o $ENV{MT_STREAM_HSW}/f_jump_ahead_coeff/gf2xe.o
        $ENV{MT_STREAM_HSW}/f_jump_ahead_coeff/f_get_coeff.o)
endif()

target_link_libraries(stochastic-mhd.exec PRIVATE stochastic MPI::MPI_Fortran flap)

add_executable(organize_mhd_data.exec organize_mhd_data.f90)
target_link_libraries(organize_mhd_data.exec PRIVATE stochastic MPI::MPI_Fortran flap)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install (TARGETS stochastic-mhd.exec DESTINATION ${CMAKE_INSTALL_PREFIX})

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install (TARGETS organize_mhd_data.exec DESTINATION ${CMAKE_INSTALL_PREFIX})

#include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# # For a large number of source files you can create it in a simpler way
# # using file() function:
# file(GLOB analysis_SOURCES *.f90)

# Create an executable files.
add_executable(stochastic-mhd.exec stochastic-mhd.f90
    $ENV{MT_STREAM}/mt_stream.o $ENV{MT_STREAM}/f_jump_ahead_coeff/gf2xe.o
    $ENV{MT_STREAM}/f_jump_ahead_coeff/f_get_coeff.o)

target_link_libraries(stochastic-mhd.exec PRIVATE stochastic MPI::MPI_Fortran flap)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install (TARGETS stochastic-mhd.exec DESTINATION ${CMAKE_INSTALL_PREFIX})

#include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# # For a large number of source files you can create it in a simpler way
# # using file() function:
# file(GLOB analysis_SOURCES *.f90)

# Create an executable files.
add_executable(stochastic-2dmhd stochastic-2dmhd.f90)
target_link_libraries(stochastic-2dmhd ${MPI_Fortran_LIBRARIES} stochastic)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install (TARGETS stochastic-2dmhd DESTINATION ${CMAKE_INSTALL_PREFIX})

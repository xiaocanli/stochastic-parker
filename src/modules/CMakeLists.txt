#include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add libraries.
add_library(stochastic STATIC constants.f90 mpi_module.f90 read_config.f90
    mhd_data_sli.f90 particle_module.f90)
#target_link_libraries(stochastic)

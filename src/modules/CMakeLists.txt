#include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add libraries.

add_library(hdf5_io hdf5_io.f90 constants.f90 mpi_module.f90)
target_link_libraries(hdf5_io ${HDF5_Fortran_LIBRARIES})

add_library(stochastic STATIC constants.f90 mpi_module.f90 read_config.f90
    mhd_data_sli.f90 particle_module.f90 random_number_generator.f90
    simulation_setup.f90 mpi_io.f90 mhd_data_parallel.f90 mhd_config.f90)
target_link_libraries(stochastic hdf5_io)

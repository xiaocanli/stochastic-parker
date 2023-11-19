cmake_minimum_required(VERSION 3.9)

project("stochastic_parker")

# omit environment variable PATH from CMake find_* use
set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH false)

#------------------------------------------------------------------------------#
# Quiet newer versions of cmake about HDF5_ROOT
#------------------------------------------------------------------------------#
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

enable_language (Fortran)

#~----------------------------------------------------------------------------~#
# MPI libraries
#~----------------------------------------------------------------------------~#

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

#------------------------------------------------------------------------------#
# Add build options
#------------------------------------------------------------------------------#

option(USE_OPENMP "Use OpenMP" OFF)
option(USE_AVX512 "Use AVX512" OFF)

#------------------------------------------------------------------------------#
# Fortran compiler and flags
#------------------------------------------------------------------------------#

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -march=native -fimplicit-none \
    -Wline-truncation -fwhole-file -std=gnu -fall-intrinsics -fPIC -cpp \
    -fallow-argument-mismatch")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fimplicit-none -Wall \
    -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing \
    -Wunused-parameter -fwhole-file -fcheck=all -std=gnu -pedantic \
    -fbacktrace -fPIC -cpp -fallow-argument-mismatch")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-xHost -O3 -g -no-prec-div -fPIC \
    -cpp -qopt-report=5 -qopt-report-phase=vec")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all \
    -check bounds -check uninit -debug all -check arg_temp_created \
    -fpe0 -ftrapuv -fPIC -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "ftn.*")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
      set (CMAKE_Fortran_FLAGS_RELEASE "-g -no-prec-div -fPIC -cpp")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -march=native -fimplicit-none \
        -Wline-truncation -fwhole-file -std=gnu -fall-intrinsics \
        -fallow-argument-mismatch -fPIC -cpp")
  endif()
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fPIC -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
endif()

#------------------------------------------------------------------------------#
# Add options for building with a threading model.
#------------------------------------------------------------------------------#

if(USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DUSE_OPENMP)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

if(USE_AVX512)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xMIC-AVX512")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -xCORE-AVX2")
endif(USE_AVX512)

#------------------------------------------------------------------------------#
# HDF5 library
#------------------------------------------------------------------------------#

set(HDF5_USE_STATIC_LIBRARIES OFF)
find_package(HDF5 COMPONENTS Fortran REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_Fortran_INCLUDE_DIRS})
link_directories(${HDF5_Fortran_LIBRARIES})
MESSAGE(STATUS ${HDF5_Fortran_LIBRARIES})

#~----------------------------------------------------------------------------~#
# FLAP to deal with command line arguments
#~----------------------------------------------------------------------------~#

include(cmake/submodules.cmake)
set(FLAP_DIR ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/FLAP)
include_directories(${FLAP_DIR}/static/mod)
link_directories(${FLAP_DIR}/static)

# use mt_stream to generate random numbers
include_directories($ENV{MT_STREAM})

# Make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(LIB ${CMAKE_SOURCE_DIR}/lib)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

# Have the .mod files placed in the lib folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# Add subdirectories
add_subdirectory(src/modules ${LIB})
add_subdirectory(src/programs ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)

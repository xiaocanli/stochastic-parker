if(EXISTS "${PROJECT_SOURCE_DIR}/.gitmodules")
message(STATUS "Updating submodules to their latest/fixed versions")
message(STATUS "(this can take a while, please be patient)")

### set the direcory where the submodules live
set(GIT_SUBMODULES_DIRECTORY src/third_party)

### set the directory names of the submodules
set(GIT_SUBMODULES FLAP)

### set each submodules's commit or tag that is to be checked out
### (leave empty if you want master)
# set(GIT_SUBMODULE_VERSION_serialization v0.9.1)
# set(GIT_SUBMODULE_VERSION_testing       4f8c51c)
# set(GIT_SUBMODULE_VERSION_make_unique   1.0.0)

### First, get all submodules in
if(${GIT_SUBMODULES_CHECKOUT_QUIET})
    execute_process(
        COMMAND             git submodule update --init --recursive
        WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}
        OUTPUT_QUIET
        ERROR_QUIET
    )
else()
    execute_process(
        COMMAND             git submodule update --init --recursive
        WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}
    )
endif()

### Then, checkout each submodule to the specified commit
# Note: Execute separate processes here, to make sure each one is run,
# should one crash (because of branch not existing, this, that ... whatever)
foreach(GIT_SUBMODULE ${GIT_SUBMODULES})

    if( "${GIT_SUBMODULE_VERSION_${GIT_SUBMODULE}}" STREQUAL "" )
        message(STATUS "no specific version given for submodule ${GIT_SUBMODULE}, checking out master")
        set(GIT_SUBMODULE_VERSION_${GIT_SUBMODULE} "master")
    endif()

    if(${GIT_SUBMODULES_CHECKOUT_QUIET})
        execute_process(
            COMMAND             git checkout ${GIT_SUBMODULE_VERSION_${GIT_SUBMODULE}}
            WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/${GIT_SUBMODULE}
            OUTPUT_QUIET
            ERROR_QUIET
        )
    else()
        message(STATUS "checking out ${GIT_SUBMODULE}'s commit/tag ${GIT_SUBMODULE_VERSION_${GIT_SUBMODULE}}")
        execute_process(
            COMMAND             git checkout ${GIT_SUBMODULE_VERSION_${GIT_SUBMODULE}}
            WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/${GIT_SUBMODULE}
        )
    endif()

endforeach(${GIT_SUBMODULE})


foreach(GIT_SUBMODULE ${GIT_SUBMODULES})
    if( ${GIT_SUBMODULE} STREQUAL "FLAP" )
        if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
            execute_process(
                COMMAND             FoBiS.py build -mode static-gnu
                WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/${GIT_SUBMODULE}
            )
        elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
            execute_process(
                COMMAND             FoBiS.py build -mode static-intel
                WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/${GIT_SUBMODULE}
            )
        elseif (Fortran_COMPILER_NAME MATCHES "ftn.*")
            if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
                execute_process(
                    COMMAND             FoBiS.py build -mode static-intel
                    WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/${GIT_SUBMODULE}
                )
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                execute_process(
                    COMMAND             FoBiS.py build -mode static-gnu
                    WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/${GIT_SUBMODULES_DIRECTORY}/${GIT_SUBMODULE}
                )
            endif()
        endif()
    endif()
endforeach(${GIT_SUBMODULE})

endif()

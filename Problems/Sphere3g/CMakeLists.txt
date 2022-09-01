cmake_minimum_required(VERSION 3.12)

#-------------------- Dependencies
if (NOT DEFINED PDES_DIR)
    if (NOT (DEFINED ENV{PDES_DIR}))
        message(FATAL_ERROR "***** PDES_DIR is not set *****")
    else()
        set(PDES_DIR "$ENV{PDES_DIR}")
    endif()
endif()
message(STATUS "PDES_DIR set to ${PDES_DIR}")

#-------------------- Define project
set(TARGET Sphere3g)
project(${TARGET} C CXX)

#-------------------- Set CMake variables
if (NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
endif()

#-------------------- Include macros
include("${PDES_DIR}/Macros/downstream.cmake")

#-------------------- Create executable
add_executable(${TARGET} "run.cc")
target_link_libraries(${TARGET} ${PDELibs})

#-------------------- Write Makefile to root directory
set(CMAKE_VERBOSE_MAKEFILE OFF)
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile
        "subsystem:\n" "\t$(MAKE) -C build \n\n"
        "clean:\n\t$(MAKE) -C build clean\n")

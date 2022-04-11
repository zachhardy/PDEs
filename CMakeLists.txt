cmake_minimum_required(VERSION 3.12)

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#-------------------- Dependencies
if (NOT DEFINED PETSC_ROOT)
    if (NOT (DEFINED ENV{PETSC_ROOT}))
        message(FATAL_ERROR "***** PETSC_ROOT is not set *****")
    else()
        set(PETSC_ROOT "$ENV{PETSC_ROOT}")
    endif()
endif()
message(STATUS "PETSC_ROOT set to ${PETSC_ROOT}")


if (NOT DEFINED EIGEN_DIR)
    if (NOT (DEFINED ENV{EIGEN_DIR}))
        message(FATAL_ERROR "***** EIGEN_DIR is not set *****")
    else()
        set(EIGEN_DIR "$ENV{EIGEN_DIR}" )
    endif()
endif()
message(STATUS "EIGEN_DIR set to ${EIGEN_DIR}")

#-------------------- Define project
set(TARGET PDEs)
project(${TARGET} C CXX)
set(MODULES "Modules")

#-------------------- CMake variables
if (NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
    set(CONFIG_FILE_NAME "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.cmake")
    file(WRITE "${CONFIG_FILE_NAME}"
            "# Environment configuration file for PDEs.\n")
endif()

#-------------------- Include macros
include_directories(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

#-------------------- Find packages
find_package(MPI)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

#-------------------- Include directories
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
include_directories(SYSTEM "${PETSC_ROOT}/include")
include_directories(SYSTEM "${EIGEN_DIR}/include/eigen3")

include_directories("${PROJECT_SOURCE_DIR}/PDEs")
include_directories("${PROJECT_SOURCE_DIR}/PDEs/Grid")
include_directories("${PROJECT_SOURCE_DIR}/PDEs/SpatialDiscretization")
include_directories("${PROJECT_SOURCE_DIR}/PDEs/Material")

include_directories("${PROJECT_SOURCE_DIR}/Modules")


#-------------------- Source directories
set(SOURCES)
add_subdirectory("${PROJECT_SOURCE_DIR}/PDEs/Grid")
add_subdirectory("${PROJECT_SOURCE_DIR}/PDEs/SpatialDiscretization")
add_subdirectory("${PROJECT_SOURCE_DIR}/PDEs/Material")

add_subdirectory("${PROJECT_SOURCE_DIR}/Modules")


#-------------------- Define targets
add_library(PDELib STATIC ${SOURCES})
target_link_libraries(PDELib ${MPI_CXX_LIBRARIES})
add_executable(${TARGET} PDEs/main.cc)
target_link_libraries(${TARGET} PDELib Eigen3::Eigen)


#-------------------- Write Makefile to root directory
set(CMAKE_VERBOSE_MAKEFILE OFF)
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile
           "subsystem:\n" "\t$(MAKE) -C build \n\n"
           "clean:\n\t$(MAKE) -C build clean\n")

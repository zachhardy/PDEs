message (STATUS "Loading Downstream.cmake")

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#-------------------- Dependencies
if (NOT DEFINED PDES_DIR)
    if (NOT (DEFINED ENV{PDES_DIR}))
        message(FATAL_ERROR "***** PDES_DIR is not set *****")
    else()
        set(PDES_DIR "$ENV{PDES_DIR}")
    endif()
endif()
message(STATUS "PDES_DIR set to ${PDES_DIR}")

include("${PDES_DIR}/bin/config.cmake")

if (NOT DEFINED PETSC_ROOT)
    if (NOT (DEFINED ENV{PETSC_ROOT}))
        message(FATAL_ERROR "***** PETSC_ROOT is not set *****")
    else()
        set(PETSC_ROOT "$ENV{PETSC_ROOT}")
    endif()
endif()
message(STATUS "PETSC_ROOT set to ${PETSC_ROOT}")

find_package(MPI)

#-------------------- Macros
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

#-------------------- Includes
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
include_directories(SYSTEM "${PETSC_ROOT}/include")

include_directories("${PDES_DIR}/PDEs")
include_directories("${PDES_DIR}/PDEs/Grid")
include_directories("${PDES_DIR}/PDEs/Math")
include_directories("${PDES_DIR}/PDEs/Physics")
include_directories("${PDES_DIR}/PDEs/Timer")

include_directories("${PDES_DIR}/Modules")

#-------------------- Libraries
link_directories("${PETSC_ROOT}/lib")
link_directories("${PDES_DIR}/lib")

set(PDELibs ${MPI_CXX_LIBRARIES} petsc)
set(PDELibs ${PDELibs} PDELib)

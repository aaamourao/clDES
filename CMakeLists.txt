# Cldes
#
# GPL v3
# Copyright (c) Adriano Madc0ww Mourão <adrianodasho@prontonmail.com>

cmake_minimum_required (VERSION 3.10)
project (clDES VERSION 0.1.0 DESCRIPTION
    "Discrete Event Systems implementation and operations")

include(ExternalProject)
set (CMAKE_CXX_STANDARD 11)

# The version number
set (clDES_VERSION_MAJOR 0)
set (clDES_VERSION_MINOR 3)
set (clDES_VERSION_PATCH 0)

# OpenCL and Boost packages are required
find_package(OpenCL 1.2 REQUIRED)
find_package(Boost 1.64 REQUIRED)

# Vars for configuring outputs dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libclDES)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libclDES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# # Vars for cloning and installing ViennaCL library
# set(viennacl_PREFIX
#     "${CMAKE_CURRENT_BINARY_DIR}/external_dep/ViennaCL-1.7.1-repo")
# set(viennacl_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/external_dep/libviennacl")
# set(viennacl_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${viennacl_INSTALL_DIR}")
#
# # Clone ViennaCL from GIT repo
# ExternalProject_Add(project_viennacl
#     GIT_REPOSITORY https://github.com/madc0ww/ViennaCL.git
#     PREFIX ${viennacl_PREFIX}
#     INSTALL_DIR ${viennacl_INSTALL_DIR}
#     CMAKE_ARGS ${viennacl_CMAKE_ARGS}
# )

# Configure a header file to pass some of the CMake settings
# to the source code
configure_file (
    "${PROJECT_SOURCE_DIR}/clDESConfig.h.in"
    "${PROJECT_BINARY_DIR}/include/clDESConfig.h"
)

# Add the binary tree to the search path for include files
# so that we will find clDESConfig.h
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include")

# Add subdirectories
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/benchmark")

enable_testing()

add_test(basics bin/tests/basics)
#add_test(ublas_basics bin/tests/ublas_basics)
add_test(sync bin/tests/sync)
add_test(sync_cpu bin/tests/sync_cpu)
add_test(direct_sync bin/tests/direct_sync)
add_test(direct_sync_cpu bin/tests/direct_sync_cpu)
add_test(supervisorsynth_cpu bin/tests/supervisorsynth_cpu)

cmake_minimum_required (VERSION 3.10)
project (clDES)

include(ExternalProject)

# The version number
set (clDES_VERSION_MAJOR 0)
set (clDES_VERSION_MINOR 3)
set (clDES_VERSION_PATCH 0)

# OpenCL package is required
find_package(OpenCL REQUIRED)

# Vars for configuring outputs dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Vars for cloning and installing ViennaCL library
set(viennacl_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/lib/ViennaCL-1.7.1")
set(viennacl_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/lib/ViennaCL-1.7.1")
set(viennacl_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${viennacl_INSTALL_DIR}")

# Clone ViennaCL from GIT repo
ExternalProject_Add(project_viennacl
    GIT_REPOSITORY https://github.com/madc0ww/ViennaCL.git
    PREFIX ${viennacl_PREFIX}
    INSTALL_DIR ${viennacl_INSTALL_DIR}
    CMAKE_ARGS ${viennacl_CMAKE_ARGS}
)

# Configure a header file to pass some of the CMake settings
# to the source code
configure_file (
    "${PROJECT_SOURCE_DIR}/include/clDESConfig.h.in"
    "${PROJECT_BINARY_DIR}/include/clDESConfig.h"
)

# Add the binary tree to the search path for include files
# so that we will find clDESConfig.h
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/lib/ViennaCL-1.7.1/include")

add_executable(clDES src/clDES.cpp)

# Link OpenCL library
target_link_libraries(clDES OpenCL::OpenCL)

# Add install targets
# install(TARGETS clDES DESTINATION "${PROJECT_SOURCE_DIR}/bin")

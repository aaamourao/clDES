# Cldes
#
# GPL v3
# Copyright (c) Adriano Madc0ww Mourão <adrianodasho@prontonmail.com>

# Library sources
set(libclDESsrcs desystem.cpp)

# Add executables
add_library(libclDES SHARED ${libclDESsrcs})

# Set library version
set_target_properties(libclDES PROPERTIES VERSION ${PROJECT_VERSION})

# Link OpenCL library
target_link_libraries(libclDES OpenCL::OpenCL)

# Add install targets
install(TARGETS libclDES DESTINATION "${PROJECT_SOURCE_DIR}/libclDES/lib")

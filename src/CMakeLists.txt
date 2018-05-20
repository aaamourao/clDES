# Cldes
#
# GPL v3
# Copyright (c) Adriano Madc0ww Mourão <adrianodasho@prontonmail.com>

# Library sources
set(libclDESSrcs
    operations.cpp
    desystem.cpp
    desystemcl.cpp
    transition_proxy.cpp
    oclbackend.cpp
)

# Library public headers
set(headersDir "${PROJECT_SOURCE_DIR}/include")

# Set local install dir
set(localInstallDir "${PROJECT_SOURCE_DIR}/libclDES/")

# Add executables
add_library(libclDES SHARED ${libclDESSrcs})

# Set library version
set_target_properties(libclDES PROPERTIES VERSION ${PROJECT_VERSION})

# Link OpenCL library
target_link_libraries(libclDES OpenCL::OpenCL)


# Add install targets
install(TARGETS libclDES
    LIBRARY DESTINATION "${localInstallDir}/lib"
)
install(DIRECTORY ${headersDir}
    DESTINATION "${localInstallDir}"
)

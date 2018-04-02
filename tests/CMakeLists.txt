# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add executables
add_executable(basics ./basics.cpp)

# Link libraries
target_link_libraries(basics libclDES OpenCL::OpenCL)

# Link OpenCL library
# target_link_libraries(viennacl OpenCL::OpenCL)

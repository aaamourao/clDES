# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

# Add executables
add_executable(benchmark ./benchmark.cpp)

# Link libraries
target_link_libraries(basics libclDES OpenCL::OpenCL)

# Link OpenCL library
# target_link_libraries(viennacl OpenCL::OpenCL)

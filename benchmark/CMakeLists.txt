# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmark)

# Add executables
add_executable(benchmark ./benchmark.cpp ./clustertools.cpp)

# Link libraries
target_link_libraries(benchmark libclDES OpenCL::OpenCL)

# Link OpenCL library
# target_link_libraries(viennacl OpenCL::OpenCL)

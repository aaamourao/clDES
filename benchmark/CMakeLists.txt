# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmark)

# Add executables
add_executable(benchmark_clustertool ./benchmark_clustertool.cpp ./clustertools.cpp)
add_executable(benchmark_fsm ./benchmark_fsm.cpp)

# Link libraries
target_link_libraries(benchmark_clustertool libclDES)# OpenCL::OpenCL)
target_link_libraries(benchmark_fsm libclDES)# OpenCL::OpenCL)

# Link OpenCL library
# target_link_libraries(viennacl OpenCL::OpenCL)

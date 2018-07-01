# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmark)

# Add executables
add_executable(benchmark_clustertool2 ./benchmark_clustertool2.cpp)
add_executable(benchmark_clustertool3 ./benchmark_clustertool3.cpp)
add_executable(benchmark_clustertool4 ./benchmark_clustertool4.cpp)
add_executable(benchmark_clustertool5 ./benchmark_clustertool5.cpp)
add_executable(benchmark_clustertool6 ./benchmark_clustertool6.cpp)
add_executable(benchmark_clustertool7 ./benchmark_clustertool7.cpp)
add_executable(benchmark_fsm ./benchmark_fsm.cpp)
add_executable(benchmark_lazyfsm ./benchmark_lazyfsm.cpp)
add_executable(benchmark_lazyclustertool6 ./benchmark_lazyclustertool6.cpp)
add_executable(benchmark_lazyclustertool7 ./benchmark_lazyclustertool7.cpp)

add_executable(benchmark_gpu_clustertool2 ./benchmark_gpu_clustertool2.cpp)
add_executable(benchmark_gpu_clustertool6 ./benchmark_gpu_clustertool6.cpp)
add_executable(benchmark_gpu_clustertool7 ./benchmark_gpu_clustertool7.cpp)
add_executable(benchmark_gpu_fsm ./benchmark_gpu_fsm.cpp)
# Link libraries
target_link_libraries(benchmark_gpu_clustertool2 OpenCL::OpenCL)
target_link_libraries(benchmark_gpu_clustertool6 OpenCL::OpenCL)
target_link_libraries(benchmark_gpu_clustertool7 OpenCL::OpenCL)
target_link_libraries(benchmark_gpu_fsm OpenCL::OpenCL)

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

# Link libraries
# target_link_libraries(benchmark_clustertool2 Qt5::Core)# OpenCL::OpenCL)
# target_link_libraries(benchmark_clustertool3 Qt5::Core)# OpenCL::OpenCL)
# target_link_libraries(benchmark_clustertool4 Qt5::Core)# OpenCL::OpenCL)
# target_link_libraries(benchmark_clustertool5 Qt5::Core)# OpenCL::OpenCL)
# target_link_libraries(benchmark_clustertool6 Qt5::Core)# OpenCL::OpenCL)
# target_link_libraries(benchmark_clustertool7 Qt5::Core)# OpenCL::OpenCL)
# target_link_libraries(benchmark_fsm Qt5::Core)# OpenCL::OpenCL)

# Link OpenCL library
# target_link_libraries(viennacl OpenCL::OpenCL)

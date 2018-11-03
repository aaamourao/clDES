# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

# Add executables
add_executable(basics ./basics.cpp)
# add_executable(basics_gpu ./basics_gpu.cpp)
add_executable(direct_sync_cpu ./direct_sync_cpu.cpp)
add_executable(supervisorsynth_cpu ./supervisorsynth_cpu.cpp)
add_executable(sync_supervisorsynth_cpu ./sync_supervisorsynth_cpu.cpp)
add_executable(superproxy ./superproxy.cpp)

# Link libraries
# target_link_libraries(basics_gpu OpenCL::OpenCL)

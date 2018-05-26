# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

# Add executables
add_executable(basics ./basics.cpp)
#add_executable(ublas_basics ./ublas_basics.cpp)
#add_executable(sync ./sync.cpp)
add_executable(sync_cpu ./sync_cpu.cpp)
#add_executable(direct_sync ./direct_sync.cpp)
add_executable(direct_sync_cpu ./direct_sync_cpu.cpp)
add_executable(supervisorsynth_cpu ./supervisorsynth_cpu.cpp)
add_executable(sync_supervisorsynth_cpu ./sync_supervisorsynth_cpu.cpp)

# Link libraries
target_link_libraries(basics libclDES OpenCL::OpenCL)
#target_link_libraries(ublas_basics libclDES OpenCL::OpenCL)
#target_link_libraries(sync libclDES OpenCL::OpenCL)
target_link_libraries(sync_cpu libclDES OpenCL::OpenCL)
#target_link_libraries(direct_sync libclDES OpenCL::OpenCL)
target_link_libraries(direct_sync_cpu libclDES OpenCL::OpenCL)
target_link_libraries(supervisorsynth_cpu libclDES OpenCL::OpenCL)
target_link_libraries(sync_supervisorsynth_cpu libclDES OpenCL::OpenCL)

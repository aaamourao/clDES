# Set output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

# Add executables
add_executable(basics ./basics.cpp)
add_executable(ublas_basics ./ublas_basics.cpp)
add_executable(sync ./sync.cpp)
add_executable(direct_sync ./direct_sync.cpp)
#add_executable(matrop ./matrop.cpp)

# Link libraries
target_link_libraries(basics libclDES OpenCL::OpenCL)
target_link_libraries(ublas_basics libclDES OpenCL::OpenCL)
target_link_libraries(sync libclDES OpenCL::OpenCL)
target_link_libraries(direct_sync libclDES OpenCL::OpenCL)
#target_link_libraries(matrop libclDES OpenCL::OpenCL)

# Link OpenCL library
# target_link_libraries(viennacl OpenCL::OpenCL)

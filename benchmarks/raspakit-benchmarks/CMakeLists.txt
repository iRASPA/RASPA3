project(benchmarks_raspakit)

add_executable(${PROJECT_NAME}
  interactions.cpp
  pow.cpp
  main.cpp)

if(CONDA)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC
    libs::raspakit
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
    OpenCL::OpenCL
    ZLIB::ZLIB
    -lpthread
    benchmark::benchmark)
elseif(LINUX)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC
    libs::raspakit
    -static-libgcc
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
    OpenCL::OpenCL
    -nostdlib++
    -l:libc++.a
    -l:libc++abi.a
    -lm
    -ldl
    ${ZLIB_LIBRARIES}
    -lpthread
    benchmark::benchmark)
elseif(WIN32)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC
    libs::raspakit
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    OpenMP::OpenMP_CXX
    OpenCL::OpenCL
    -nostdlib++
    ${ZLIB_LIBRARIES}
    benchmark::benchmark)
elseif(APPLE)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC
    libs::raspakit
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${TORCH_LIBRARIES}
    $ENV{HOMEBREW_PREFIX}/lib/libhdf5_cpp.a
    $ENV{HOMEBREW_PREFIX}/lib/libhdf5.a
    $ENV{HOMEBREW_PREFIX}/lib/libaec.a
    $ENV{HOMEBREW_PREFIX}/lib/libsz.a
    -nostdlib++
    $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++experimental.a
    $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++.a
    $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++abi.a
    $ENV{HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
    OpenCL::OpenCL
    ${ZLIB_LIBRARIES}
    -lpthread
    benchmark::benchmark)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/benchmarks COMPONENT ${PROJECT_NAME})

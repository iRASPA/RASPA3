add_executable(benchmarks_raspakit
               benchmark1.cpp
               benchmark2.cpp
               main.cpp)


if (LINUX)
  target_link_libraries(benchmarks_raspakit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        benchmark::benchmark
                        -nostdlib++
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -ldl
                        -lz
                        -lpthread
)
elseif (Win32)
  target_link_libraries(benchmarks_raspakit
                        PUBLIC
                        libs::raspakit
                        -static-libstdc++
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -lpthread
)
elseif (APPLE)
  target_link_libraries(benchmarks_raspakit PRIVATE
    libs::raspakit
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    benchmark::benchmark
    -nostdlib++
    ${HOMEBREW_PREFIX}/opt/llvm//lib/c++/libc++experimental.a
    ${HOMEBREW_PREFIX}/opt/llvm//lib/c++/libc++.a
    ${HOMEBREW_PREFIX}/opt/llvm//lib/c++/libc++abi.a
    ${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
    -lpthread
  )
endif()

install(TARGETS benchmarks_raspakit DESTINATION share/raspa3/benchmarks COMPONENT benchmarks_raspakit)

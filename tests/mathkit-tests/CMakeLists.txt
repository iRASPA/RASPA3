add_executable(unit_tests_mathkit
               double3x3.cpp
               main.cpp)

if (CONDA)
  target_link_libraries(unit_tests_mathkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        GTest::gtest
                        GTest::gtest_main
                        ZLIB::ZLIB
                        -lpthread
)
elseif (LINUX)
  target_link_libraries(unit_tests_mathkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        GTest::gtest 
                        GTest::gtest_main 
                        -nostdlib++
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -ldl
                        ${ZLIB_LIBRARIES}
                        -lpthread
)
elseif (WIN32 AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "AMD64")
  target_link_libraries(unit_tests_mathkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        GTest::gtest 
                        GTest::gtest_main 
                        ${llvm_libs}
                        -nostdlib++
)
elseif (WIN32 AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "ARM64")
  target_link_libraries(unit_tests_mathkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        C:/vcpkg/installed/arm64-windows/lib/libf2c.lib
                        ${HDF5_CXX_LIBRARIES}
                        GTest::gtest 
                        GTest::gtest_main 
                        ${llvm_libs}
                        -nostdlib++
)
elseif (APPLE)
  target_link_libraries(unit_tests_mathkit PRIVATE
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        $ENV{HOMEBREW_PREFIX}/lib/libhdf5_cpp.a
                        $ENV{HOMEBREW_PREFIX}/lib/libhdf5.a
                        $ENV{HOMEBREW_PREFIX}/lib/libaec.a
                        $ENV{HOMEBREW_PREFIX}/lib/libsz.a
                        GTest::gtest
                        GTest::gtest_main
                        -nostdlib++
                        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++experimental.a
                        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++.a
                        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++abi.a
                        $ENV{HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
                        ${ZLIB_LIBRARIES}
                        -lpthread
  ) 
endif()

gtest_add_tests(TARGET unit_tests_mathkit TEST_SUFFIX .noArgs TEST_LIST   noArgsTests)

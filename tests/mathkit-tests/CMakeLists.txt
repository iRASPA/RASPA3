add_executable(unit_tests_mathkit
               double3x3.cpp
               main.cpp)

if (LINUX)
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
                        GTest::gmock 
                        GTest::gmock_main 
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
                        ${BREW_HDF5_DIR}/lib/libhdf5_cpp.a
                        ${BREW_HDF5_DIR}/lib/libhdf5.a
                        ${BREW_LIBAEC_DIR}/lib/libaec.a
                        ${BREW_LIBAEC_DIR}/lib/libsz.a
                        GTest::gtest
                        GTest::gtest_main
                        GTest::gmock
                        GTest::gmock_main
                        -nostdlib++
                        ${BREW_LLVM_DIR}/lib/c++/libc++experimental.a
                        ${BREW_LLVM_DIR}/lib/c++/libc++.a
                        ${BREW_LLVM_DIR}/lib/c++/libc++abi.a
                        ${BREW_LIBOMP_DIR}/lib/libomp.a
                        ${ZLIB_LIBRARIES}
                        -lpthread
  ) 
endif()

gtest_discover_tests(unit_tests_mathkit)

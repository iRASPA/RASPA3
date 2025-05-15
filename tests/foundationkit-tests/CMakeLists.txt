add_executable(unit_tests_foundationkit
               archive.cpp 
               main.cpp)


if (CONDA)
  target_link_libraries(unit_tests_foundationkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        voro++
                        ${TORCH_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        OpenCL::OpenCL
                        GTest::gtest
                        GTest::gtest_main
                        ZLIB::ZLIB
                        -lpthread)
elseif (LINUX)
  target_link_libraries(unit_tests_foundationkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        voro++
                        ${TORCH_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        OpenCL::OpenCL
                        GTest::gtest
                        GTest::gtest_main
                        -nostdlib++
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -ldl
                        ${ZLIB_LIBRARIES}
                        -lpthread)
elseif (WIN32)
  target_link_libraries(unit_tests_foundationkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        voro++
                        ${TORCH_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        OpenCL::OpenCL
                        GTest::gtest
                        GTest::gtest_main
                        -nostdlib++
                        ${ZLIB_LIBRARIES})
elseif (APPLE)
  target_link_libraries(unit_tests_foundationkit PRIVATE
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        $ENV{HOMEBREW_PREFIX}/lib/libhdf5_cpp.a
                        $ENV{HOMEBREW_PREFIX}/lib/libhdf5.a
                        $ENV{HOMEBREW_PREFIX}/lib/libaec.a
                        $ENV{HOMEBREW_PREFIX}/lib/libsz.a
                        voro++
                        ${TORCH_LIBRARIES}
                        GTest::gtest
                        GTest::gtest_main
                        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++experimental.a
                        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++.a
                        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++abi.a
                        $ENV{HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
                        OpenCL::OpenCL
                        ${ZLIB_LIBRARIES}
                        -lpthread)
endif()

gtest_add_tests(TARGET unit_tests_foundationkit TEST_SUFFIX .noArgs TEST_LIST   noArgsTests)

if(BUILD_MAC_PACKAGE)   
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"
      \"tests/foundationkit-tests/unit_tests_foundationkit\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_foundationkit DESTINATION share/raspa3/tests COMPONENT unit_tests)


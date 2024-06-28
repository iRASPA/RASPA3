add_executable(unit_tests_foundationkit
               archive.cpp 
               main.cpp)


if (LINUX)
  target_link_libraries(unit_tests_foundationkit
                        PUBLIC
                        libs::raspakit
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
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
                        -lpthread
)
elseif (Win32)
  target_link_libraries(unit_tests_foundationkit
                        PUBLIC
                        libs::raspakit
                        -static-libstdc++
                        ${BLAS_LIBRARIES}
                        ${LAPACK_LIBRARIES}
                        OpenMP::OpenMP_CXX
                        -l:libc++.a
                        -l:libc++abi.a
                        -lm
                        -lpthread
)
elseif (APPLE)
target_link_libraries(unit_tests_foundationkit PRIVATE
                      libs::raspakit
                      ${BLAS_LIBRARIES}
                      ${LAPACK_LIBRARIES}
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main
                      -nostdlib++
                      ${LLVM_LIBS}/c++/libc++experimental.a
                      ${LLVM_LIBS}/c++/libc++experimental.a
                      ${LLVM_LIBS}/c++/libc++.a
                      ${LLVM_LIBS}/c++/libc++abi.a
                      ${LLVM_LIBS}/libunwind.a
                      ${LLVM_LIBS}/../../../libomp/${LLVM_VERSION}/lib/libomp.a
                      -lpthread
)
endif()

gtest_discover_tests(unit_tests_foundationkit)

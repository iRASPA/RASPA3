add_executable(unit_tests_raspakit
  cbmc_interactions.cpp
  energy_decomposition.cpp
  gradients.cpp
  integrators.cpp
  pbc.cpp
  spacegroup.cpp
  dudlambda.cpp
  ewald.cpp
  insertion_deletion.cpp
  mc.cpp
  pressure.cpp
  static_energy.cpp
  rigid_gradient.cpp
  electrostatic_potential.cpp
  electrostatic_field.cpp
  electrostatic_polarization.cpp
  nist_water_1.cpp
  nist_water_2.cpp
  nist_water_3.cpp
  nist_water_4.cpp
  hessian_inter_lennard_jones.cpp
  hessian_inter_real_ewald.cpp
  tricubic_derivative_inter_lennard_jones.cpp
  triquintic_derivative_inter_lennard_jones.cpp
  tricubic_derivative_inter_real_ewald.cpp
  triquintic_derivative_inter_real_ewald.cpp
  #  grids.cpp
  main.cpp)


if (CONDA)
  target_link_libraries(unit_tests_raspakit
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
  target_link_libraries(unit_tests_raspakit
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
  target_link_libraries(unit_tests_raspakit
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
  target_link_libraries(unit_tests_raspakit PRIVATE
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

gtest_add_tests(TARGET unit_tests_raspakit TEST_SUFFIX .noArgs TEST_LIST noArgsTests)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
      \"tests/raspakit-tests/unit_tests_raspakit\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_raspakit DESTINATION share/raspa3/tests COMPONENT unit_tests)


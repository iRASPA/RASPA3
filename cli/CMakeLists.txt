project(raspa3-cli)

add_executable(${PROJECT_NAME}
         main.cpp
)

target_sources(${PROJECT_NAME} PUBLIC commandline.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES commandline.ixx)

if (CONDA)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_LIBRARIES}
        voro++
        ${TORCH_LIBRARIES}
        OpenMP::OpenMP_CXX
        OpenCL::OpenCL
        ZLIB::ZLIB
        -lpthread)
elseif (LINUX)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libgcc
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        ${HDF5_LIBRARIES}
        voro++
        ${TORCH_LIBRARIES}
        OpenMP::OpenMP_CXX
        OpenCL::OpenCL
        -nostdlib++
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -ldl
        ${ZLIB_LIBRARIES}
        -lpthread)
elseif (WIN32)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        voro++
        ${TORCH_LIBRARIES}
        OpenMP::OpenMP_CXX
        OpenCL::OpenCL
        -nostdlib++
        ${ZLIB_LIBRARIES})
elseif (APPLE)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        voro++
        ${TORCH_LIBRARIES}
        $ENV{HOMEBREW_PREFIX}/lib/libhdf5_cpp.a
        $ENV{HOMEBREW_PREFIX}/lib/libhdf5.a
        $ENV{HOMEBREW_PREFIX}/lib/libaec.a 
        $ENV{HOMEBREW_PREFIX}/lib/libsz.a 
        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++experimental.a
        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++.a
        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++abi.a
        $ENV{HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
        OpenCL::OpenCL
        ${ZLIB_LIBRARIES}
        -lpthread)
endif()

if(BUILD_MAC_PACKAGE)
install(CODE "execute_process(COMMAND
    codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
    \"cli/${PROJECT_NAME}\")" COMPONENT cli)
endif()
install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT cli)

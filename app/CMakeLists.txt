project(raspa3)

add_executable(${PROJECT_NAME}
        main.cpp
)

if (CONDA)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_LIBRARIES}
        OpenMP::OpenMP_CXX
        ZLIB::ZLIB
        -lpthread
)
elseif (LINUX)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libgcc
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        ${HDF5_LIBRARIES}
        OpenMP::OpenMP_CXX
        -nostdlib++
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -ldl
        ${ZLIB_LIBRARIES}
        -lpthread)
elseif (WIN32 AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "AMD64")
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        -nostdlib++)
elseif (WIN32 AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "ARM64")
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        C:/vcpkg/installed/arm64-windows-static/lib/libf2c.lib
        ${HDF5_CXX_LIBRARIES}
        -nostdlib++)
elseif (APPLE)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        $ENV{HOMEBREW_PREFIX}/lib/libhdf5_cpp.a
        $ENV{HOMEBREW_PREFIX}/lib/libhdf5.a
        $ENV{HOMEBREW_PREFIX}/lib/libaec.a 
        $ENV{HOMEBREW_PREFIX}/lib/libsz.a 
        -nostdlib++
        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++experimental.a
        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++.a
        $ENV{HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++abi.a
        $ENV{HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
        ${ZLIB_LIBRARIES}
        -lpthread)
endif()

if(BUILD_MAC_PACKAGE)
install(CODE "execute_process(COMMAND
    codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
    \"app/${PROJECT_NAME}\")" COMPONENT app)
endif()
install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT app)

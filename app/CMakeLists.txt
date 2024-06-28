project(raspa)

add_executable(${PROJECT_NAME}
        main.cpp
)

if (LINUX)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libgcc
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        -nostdlib++
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -ldl
        -lpthread
)
elseif (Win32)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libstdc++
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        -nostdlib++
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -lpthread
)
elseif (APPLE)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        -nostdlib++
        ${LLVM_LIBS}/c++/libc++experimental.a
        ${LLVM_LIBS}/c++/libc++.a
        ${LLVM_LIBS}/c++/libc++abi.a
        ${LLVM_LIBS}/libunwind.a
        ${LLVM_LIBS}/../../../libomp/${LLVM_VERSION}/lib/libomp.a
        -lpthread
)
endif()


install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT app)

project(raspa)

add_executable(${PROJECT_NAME}
        main.cpp
)

if (LINUX)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libstdc++
        -static-libgcc
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -lpthread
)
elseif (Win32)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libstdc++
        -static-libgcc
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -lpthread
)
elseif (APPLE)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -lpthread
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
)
endif()


install(TARGETS ${PROJECT_NAME} DESTINATION bin)

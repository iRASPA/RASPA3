project(raspa)

add_executable(${PROJECT_NAME}
        main.cpp
)

target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libstdc++
        -static-libgcc
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        ${HDF5_CXX_LIBRARIES}
        -l:libc++.a
        -l:libc++abi.a
        -lm
        -lpthread
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

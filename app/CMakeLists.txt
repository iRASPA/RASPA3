project(raspa)

add_executable(${PROJECT_NAME}
	main.cpp
)

target_link_libraries(${PROJECT_NAME}
        PRIVATE
        libs::raspakit
        -lc++
        -lc++abi
        -lpthread
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        ${HDF5_CXX_LIBRARIES}
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

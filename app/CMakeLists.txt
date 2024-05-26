project(raspa)

add_executable(${PROJECT_NAME}
	main.cpp
)

target_link_libraries(${PROJECT_NAME}
        PRIVATE
        libs::raspakit
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_CXX
        -lc++
        -lc++abi
        -lpthread
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

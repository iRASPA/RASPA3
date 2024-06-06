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
        -lpthread
)

target_link_libraries(${PROJECT_NAME} PUBLIC -static-libstdc++ -l:libc++.a -l:libc++abi.a -static-libgcc)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

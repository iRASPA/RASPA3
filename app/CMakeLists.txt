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
                pybind11::module 
                ${BLAS_LIBRARIES} 
                ${LAPACK_LIBRARIES}
                Python::Python
                OpenMP::OpenMP_CXX
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

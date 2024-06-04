project(_raspa)

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(${PROJECT_NAME} MODULE pybindings.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES pybindings.ixx)
target_link_libraries(${PROJECT_NAME} PUBLIC libs::raspakit ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} OpenMP::OpenMP_CXX -lc++ -lc++abi -lpthread)

# add_library(${PROJECT_NAME} MODULE pybindings.cpp)
# target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES pybindings.ixx)
# target_link_libraries(${PROJECT_NAME} pybind11::module libs::raspakit)
# set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET "hidden")
target_compile_options(${PROJECT_NAME} PUBLIC -stdlib=libc++ $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

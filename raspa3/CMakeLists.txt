project(_raspa)
find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(${PROJECT_NAME} pybindings.cpp)

# add_library(${PROJECT_NAME} SHARED libs::raspakit)
target_link_libraries(${PROJECT_NAME} PRIVATE libs::raspakit)
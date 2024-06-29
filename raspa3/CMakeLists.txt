project(raspalib)

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(${PROJECT_NAME} MODULE)
target_sources(${PROJECT_NAME} PUBLIC pybindings.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES pybindings.ixx)

if (LINUX)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        pybind11::module
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
        pybind11::module
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
        pybind11::module
        libs::raspakit
        -nostdlib++
        ${LLVM_LIBS}/c++/libc++experimental.a
        ${LLVM_LIBS}/c++/libc++.a
        ${LLVM_LIBS}/c++/libc++abi.a
        ${LLVM_LIBS}/libunwind.a
        ${LLVM_LIBS}/../../../libomp/${LLVM_VERSION}/lib/libomp.a
        -lpthread
)
endif()

target_compile_options(${PROJECT_NAME} PUBLIC -stdlib=libc++ $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

pybind11_extension(${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION share/raspa3/lib COMPONENT python)

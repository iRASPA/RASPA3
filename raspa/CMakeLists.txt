project(raspalib)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.10.4
)

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

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
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
    -nostdlib++
    -l:libc++.a
    -l:libc++abi.a
    -lm
    -ldl
    -lpthread
  )
elseif (WIN32)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC
    libs::raspakit
    pybind11::module
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
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
    pybind11::module
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BREW_HDF5_DIR}/lib/libhdf5_cpp.a
    ${BREW_HDF5_DIR}/lib/libhdf5.a
    ${BREW_LIBAEC_DIR}/lib/libaec.a
    ${BREW_LIBAEC_DIR}/lib/libsz.a
    -nostdlib++
    ${BREW_LLVM_DIR}/lib/c++/libc++experimental.a
    ${BREW_LLVM_DIR}/lib/c++/libc++.a
    ${BREW_LLVM_DIR}/lib/c++/libc++abi.a
    ${BREW_LIBOMP_DIR}/lib/libomp.a
    ${ZLIB_LIBRARIES}
    -lpthread
  )
endif()

target_compile_options(${PROJECT_NAME} PUBLIC -stdlib=libc++ $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

if (LINUX)
  install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/lib)
elseif (APPLE)
  install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/lib)
endif()

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
    FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(${PROJECT_NAME} MODULE)
target_sources(${PROJECT_NAME} PUBLIC raspa_pybindings.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES raspa_pybindings.ixx)

if (CONDA)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        pybind11::module
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HDF5_LIBRARIES}
        OpenMP::OpenMP_CXX
        ZLIB::ZLIB
        -lpthread)
elseif (LINUX)
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
        -lpthread)
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
        -nostdlib++)
elseif (APPLE)
  target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        pybind11::module
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${HOMEBREW_PREFIX}/lib/libhdf5_cpp.a
        ${HOMEBREW_PREFIX}/lib/libhdf5.a
        ${HOMEBREW_PREFIX}/lib/libaec.a
        ${HOMEBREW_PREFIX}/lib/libsz.a
        -nostdlib++
        ${HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++experimental.a
        ${HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++.a
        ${HOMEBREW_PREFIX}/opt/llvm/lib/c++/libc++abi.a
        ${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a
        ${ZLIB_LIBRARIES}
        -lpthread)
endif()

target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)


if (CONDA)
  target_link_options(${PROJECT_NAME} PRIVATE -fuse-ld=lld -L$ENV{CONDA_PREFIX}/lib)
  execute_process ( COMMAND python -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS ${PROJECT_NAME} DESTINATION ${PYTHON_SITE_PACKAGES})
elseif (LINUX)
  install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/lib)
elseif (APPLE)
  install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/lib)
endif()

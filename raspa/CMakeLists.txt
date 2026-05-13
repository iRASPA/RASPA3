project(raspalib)

if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
  find_package(Python3 REQUIRED COMPONENTS Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.SABIModule)
endif()

if(NOT Python_FOUND)
    set(Python_FOUND ${Python3_FOUND})
    set(Python_Development_FOUND ${Python3_Development.Module_FOUND})
    set(Python_Interpreter_FOUND ${Python3_Interpreter_FOUND})
endif()

include(FetchContent)

set(NB_TEST OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind
    GIT_TAG        v2.12.0
    GIT_SHALLOW    YES
)

FetchContent_MakeAvailable(nanobind)

add_compile_definitions(Py_LIMITED_API=0x030c0000)

nanobind_build_library(nanobind-static-abi3)

add_library(${PROJECT_NAME} MODULE)
#nanobind_add_module(${PROJECT_NAME} NB_STATIC STABLE_ABI LTO raspa_pybindings.cpp)
target_sources(${PROJECT_NAME} PUBLIC raspa_pybindings.cpp)
set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_MODULE_STD ON)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_SCAN_FOR_MODULES ON)
raspa_target_link_libraries(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PRIVATE nanobind-static-abi3 Python3::SABIModule)
target_include_directories(${PROJECT_NAME} PRIVATE
  ${nanobind_SOURCE_DIR}/include
  ${nanobind_SOURCE_DIR}/ext/robin-map/include
  ${Python3_SABI_INCLUDE_DIRS}
)
include_directories(${Python3_INCLUDE_DIRS})

target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

set_target_properties(${PROJECT_NAME} PROPERTIES
  PREFIX ""
  SUFFIX "${NB_SUFFIX_S}" # Leverages nanobind's resolved extension macro
)

#if(APPLE)
#   target_link_options(${PROJECT_NAME} PRIVATE "-Wl,-undefined,dynamic_lookup")
#endif()

## on intel, also create avx2 and avx512 versions
#if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")
#
#  pybind11_add_module(${PROJECT_NAME}_avx2 MODULE)
#  target_sources(${PROJECT_NAME}_avx2 PUBLIC raspa_pybindings.cpp)
#  target_sources(${PROJECT_NAME}_avx2 PUBLIC FILE_SET CXX_MODULES FILES raspa_pybindings.ixx)
#  if(MSVC)
#    target_compile_options(${PROJECT_NAME}_avx2 PRIVATE /arch:AVX2 -DMODULE_NAME=raspalib_avx2)
#  else()
#    target_compile_options(${PROJECT_NAME}_avx2 PRIVATE -march=core-avx2 -mfma -DMODULE_NAME=raspalib_avx2)
#  endif()
#  target_compile_options(${PROJECT_NAME}_avx2 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
#  set_target_properties(${PROJECT_NAME}_avx2 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#  raspa_target_link_libraries(${PROJECT_NAME}_avx2)
#  
#  pybind11_add_module(${PROJECT_NAME}_avx512 MODULE)
#  target_sources(${PROJECT_NAME}_avx512 PUBLIC raspa_pybindings.cpp)
#  target_sources(${PROJECT_NAME}_avx512 PUBLIC FILE_SET CXX_MODULES FILES raspa_pybindings.ixx)
#  if(MSVC)
#    target_compile_options(${PROJECT_NAME}_avx512 PRIVATE /arch:AVX512 -DMODULE_NAME=raspalib_avx512)
#  else()
#    target_compile_options(${PROJECT_NAME}_avx512 PRIVATE -march=skylake-avx512 -mfma -DMODULE_NAME=raspalib_avx512)
#  endif()
#  target_compile_options(${PROJECT_NAME}_avx512 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
#  set_target_properties(${PROJECT_NAME}_avx512 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#  raspa_target_link_libraries(${PROJECT_NAME}_avx512)
#
#endif()

if(CONDA)
  target_link_options(${PROJECT_NAME} PRIVATE -L$ENV{CONDA_PREFIX}/lib)
  execute_process(COMMAND python -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS ${PROJECT_NAME} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
elseif(LINUX)
  # base the install directory on where also python-numpy is installed
  #execute_process(COMMAND python3 -c "import numpy; import os; print(os.path.dirname(os.path.dirname(numpy.__file__)))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND python3 -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS ${PROJECT_NAME} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
elseif(APPLE)
  if(PYPI)
    install(TARGETS ${PROJECT_NAME} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    install(CODE "
            execute_process(
                COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}>\"
            )
        ")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")
      install(TARGETS ${PROJECT_NAME}_avx2 DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(CODE "
              execute_process(
                  COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}_avx2>\"
              )
          ")
      install(TARGETS ${PROJECT_NAME}_avx512 DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(CODE "
              execute_process(
                  COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}_avx512>\"
              )
          ")
    endif()
    install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
  else()
    install(TARGETS ${PROJECT_NAME} DESTINATION raspalib)
    install(CODE "
            execute_process(
                COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}>\"
            )
        ")
    install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/lib)
  endif()
endif()

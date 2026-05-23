project(raspalib)

if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
  find_package(Python REQUIRED COMPONENTS Development)
else()
  find_package(Python REQUIRED COMPONENTS Interpreter Development.SABIModule)
endif()

include(FetchContent)

set(NB_TEST OFF CACHE BOOL "" FORCE)
if (Python_EXECUTABLE)
    set(NB_PYTHON_EXECUTABLE "${Python_EXECUTABLE}" CACHE INTERNAL "")
    # Also set standard CMake hints just in case
    set(Python_EXECUTABLE "${Python_EXECUTABLE}" CACHE INTERNAL "")
endif()

FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind
    GIT_TAG        v2.12.0
    GIT_SHALLOW    YES
)

FetchContent_MakeAvailable(nanobind)

# Force-enable C++20 Standard Library Modules after nanobind disables them
set(CMAKE_CXX_MODULE_STD ON CACHE BOOL "" FORCE)
set(CXX_SCAN_FOR_MODULES ON CACHE BOOL "" FORCE)
set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_VISIBILITY_INLINES_HIDDEN OFF)

add_compile_definitions(Py_LIMITED_API=0x030c0000)

nanobind_build_library(nanobind-static-abi3)

add_library(${PROJECT_NAME}_base MODULE)
target_sources(${PROJECT_NAME}_base PUBLIC raspa_pybindings.cpp)
target_compile_options(${PROJECT_NAME}_base PRIVATE -DMODULE_NAME=raspalib_base)
set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
set_property(TARGET ${PROJECT_NAME}_base PROPERTY CXX_MODULE_STD ON)
set_property(TARGET ${PROJECT_NAME}_base PROPERTY CXX_SCAN_FOR_MODULES ON)
raspa_target_link_libraries(${PROJECT_NAME}_base PRIVATE)
target_link_libraries(${PROJECT_NAME}_base PRIVATE raspakit_base nanobind-static-abi3 Python::SABIModule)
target_include_directories(${PROJECT_NAME}_base PRIVATE
  ${nanobind_SOURCE_DIR}/include
  ${nanobind_SOURCE_DIR}/ext/robin-map/include
  ${Python_SABI_INCLUDE_DIRS}
)
include_directories(${Python_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME}_base PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME}_base PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(${PROJECT_NAME}_base PROPERTIES
  PREFIX ""
  SUFFIX "${NB_SUFFIX_S}" # Leverages nanobind's resolved extension macro
)

if(PYPI)
  # on intel, also create avx, avx2 and avx512 versions
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")
  
    add_library(${PROJECT_NAME}_avx MODULE)
    target_sources(${PROJECT_NAME}_avx PUBLIC raspa_pybindings.cpp)
    if(MSVC)
      target_compile_options(${PROJECT_NAME}_avx PRIVATE /arch:AVX -DMODULE_NAME=raspalib_avx)
    else()
      # -mno-avx512f needed for Google-Colab (fixes "illegal instruction")
      target_compile_options(${PROJECT_NAME}_avx PRIVATE -march=x86-64-v2 -mno-avx512f -mavx -DMODULE_NAME=raspalib_avx)
    endif()
    set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
    set_property(TARGET ${PROJECT_NAME}_avx PROPERTY CXX_MODULE_STD ON)
    set_property(TARGET ${PROJECT_NAME}_avx PROPERTY CXX_SCAN_FOR_MODULES ON)
    raspa_target_link_libraries(${PROJECT_NAME}_avx PRIVATE)
    target_link_libraries(${PROJECT_NAME}_avx PRIVATE raspakit_avx nanobind-static-abi3 Python::SABIModule)
    target_include_directories(${PROJECT_NAME}_avx PRIVATE
      ${nanobind_SOURCE_DIR}/include
      ${nanobind_SOURCE_DIR}/ext/robin-map/include
      ${Python_SABI_INCLUDE_DIRS}
    )
    include_directories(${Python_INCLUDE_DIRS})
    target_compile_options(${PROJECT_NAME}_avx PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
    set_target_properties(${PROJECT_NAME}_avx PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_target_properties(${PROJECT_NAME}_avx PROPERTIES
      PREFIX ""
      SUFFIX "${NB_SUFFIX_S}"
    )
  
    add_library(${PROJECT_NAME}_avx2 MODULE)
    target_sources(${PROJECT_NAME}_avx2 PUBLIC raspa_pybindings.cpp)
    if(MSVC)
      target_compile_options(${PROJECT_NAME}_avx2 PRIVATE /arch:AVX2 -DMODULE_NAME=raspalib_avx2)
    else()
      # -mno-avx512f needed for Google-Colab (fixes "illegal instruction")
      target_compile_options(${PROJECT_NAME}_avx2 PRIVATE -march=x86-64-v3 -mno-avx512f -mavx2 -mfma -DMODULE_NAME=raspalib_avx2)
    endif()
    set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
    set_property(TARGET ${PROJECT_NAME}_avx2 PROPERTY CXX_MODULE_STD ON)
    set_property(TARGET ${PROJECT_NAME}_avx2 PROPERTY CXX_SCAN_FOR_MODULES ON)
    raspa_target_link_libraries(${PROJECT_NAME}_avx2 PRIVATE)
    target_link_libraries(${PROJECT_NAME}_avx2 PRIVATE raspakit_avx2 nanobind-static-abi3 Python::SABIModule)
    target_include_directories(${PROJECT_NAME}_avx2 PRIVATE
      ${nanobind_SOURCE_DIR}/include
      ${nanobind_SOURCE_DIR}/ext/robin-map/include
      ${Python_SABI_INCLUDE_DIRS}
    )
    include_directories(${Python_INCLUDE_DIRS})
    target_compile_options(${PROJECT_NAME}_avx2 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
    set_target_properties(${PROJECT_NAME}_avx2 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_target_properties(${PROJECT_NAME}_avx2 PROPERTIES
      PREFIX ""
      SUFFIX "${NB_SUFFIX_S}"
    )
  
    add_library(${PROJECT_NAME}_avx512 MODULE)
    target_sources(${PROJECT_NAME}_avx512 PUBLIC raspa_pybindings.cpp)
    if(MSVC)
      target_compile_options(${PROJECT_NAME}_avx512 PRIVATE /arch:AVX512 -DMODULE_NAME=raspalib_avx512)
    else()
      target_compile_options(${PROJECT_NAME}_avx512 PRIVATE -march=x86-64-v4 -mavx512f -mfma -mprefer-vector-width=512 -DMODULE_NAME=raspalib_avx512)
    endif()
    set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
    set_property(TARGET ${PROJECT_NAME}_avx512 PROPERTY CXX_MODULE_STD ON)
    set_property(TARGET ${PROJECT_NAME}_avx512 PROPERTY CXX_SCAN_FOR_MODULES ON)
    raspa_target_link_libraries(${PROJECT_NAME}_avx512 PRIVATE)
    target_link_libraries(${PROJECT_NAME}_avx512 PRIVATE raspakit_avx512 nanobind-static-abi3 Python::SABIModule)
    target_include_directories(${PROJECT_NAME}_avx512 PRIVATE
      ${nanobind_SOURCE_DIR}/include
      ${nanobind_SOURCE_DIR}/ext/robin-map/include
      ${Python_SABI_INCLUDE_DIRS}
    )
    include_directories(${Python_INCLUDE_DIRS})
    target_compile_options(${PROJECT_NAME}_avx512 PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
    set_target_properties(${PROJECT_NAME}_avx512 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_target_properties(${PROJECT_NAME}_avx512 PROPERTIES
      PREFIX ""
      SUFFIX "${NB_SUFFIX_S}"
    )
  
  endif()
endif()

if(CONDA)
  target_link_options(${PROJECT_NAME}_base PRIVATE -L$ENV{CONDA_PREFIX}/lib)
  execute_process(COMMAND python -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS ${PROJECT_NAME}_base DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
elseif(WIN32)
  if(PYPI)
    install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
    install(TARGETS ${PROJECT_NAME}_base DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")
      install(TARGETS ${PROJECT_NAME}_avx DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(TARGETS ${PROJECT_NAME}_avx2 DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(TARGETS ${PROJECT_NAME}_avx512 DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    endif()
  endif()
elseif(LINUX)
  if(PYPI)
    install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
    install(TARGETS ${PROJECT_NAME}_base DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")
      install(TARGETS ${PROJECT_NAME}_avx DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(TARGETS ${PROJECT_NAME}_avx2 DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(TARGETS ${PROJECT_NAME}_avx512 DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    endif()
  else()
    execute_process(COMMAND python3 -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    install(TARGETS ${PROJECT_NAME}_base DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
    install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
endif()
elseif(APPLE)
  if(PYPI)
    install(TARGETS ${PROJECT_NAME}_base DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    install(CODE "
            execute_process(
                COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}_base>\"
            )
        ")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64)")
      install(TARGETS ${PROJECT_NAME}_avx DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
      install(CODE "
              execute_process(
                  COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}_avx>\"
              )
          ")
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
    install(TARGETS ${PROJECT_NAME}_base DESTINATION raspalib)
    install(CODE "
            execute_process(
                COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"raspalib/$<TARGET_FILE_NAME:${PROJECT_NAME}_base>\"
            )
        ")
    install(TARGETS ${PROJECT_NAME}_base DESTINATION share/raspa3/lib)
  endif()
endif()

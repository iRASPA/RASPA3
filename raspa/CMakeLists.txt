project(raspalib)

# The extension module is installed inside the raspalib package and re-exported by its __init__.py, so it
# needs a name of its own: a submodule named raspalib inside a package named raspalib cannot be imported.
set(MODULE_TARGET _raspalib)

# The bindings are written with nanobind for every kind of build: the wheels for PyPI (PYPI, with the stable ABI
# when ABI3 is set), the conda packages (CONDA) and a local build that installs into the site-packages of the
# interpreter found here.

if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
  if(CONDA)
    set(Python_FIND_STRATEGY LOCATION)
    if(NOT Python_ROOT_DIR)
      set(Python_ROOT_DIR "$ENV{PREFIX}")
    endif()
  endif()
  find_package(Python REQUIRED COMPONENTS Development)
elseif(ABI3)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.SABIModule)
else()
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
endif()

include(FetchContent)

set(NB_TEST OFF CACHE BOOL "" FORCE)
if(Python_EXECUTABLE)
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

if(ABI3)
  add_compile_definitions(Py_LIMITED_API=0x030c0000)
  set(NANOBIND_LIBRARY nanobind-static-abi3)
  set(PYTHON_MODULE_TARGET Python::SABIModule)
  set(PYTHON_MODULE_INCLUDE_DIRS ${Python_SABI_INCLUDE_DIRS})
  set(PYTHON_MODULE_SUFFIX "${NB_SUFFIX_S}")
else()
  set(NANOBIND_LIBRARY nanobind-static)
  set(PYTHON_MODULE_TARGET Python::Module)
  set(PYTHON_MODULE_INCLUDE_DIRS ${Python_INCLUDE_DIRS})
  set(PYTHON_MODULE_SUFFIX "${NB_SUFFIX}")
endif()

nanobind_build_library(${NANOBIND_LIBRARY})

add_library(${MODULE_TARGET} MODULE)
target_sources(${MODULE_TARGET} PUBLIC raspa_pybindings.cpp)
target_compile_options(${MODULE_TARGET} PRIVATE -DMODULE_NAME=_raspalib)
set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
set_property(TARGET ${MODULE_TARGET} PROPERTY CXX_MODULE_STD ON)
set_property(TARGET ${MODULE_TARGET} PROPERTY CXX_SCAN_FOR_MODULES ON)
raspa_target_link_libraries(${MODULE_TARGET} PRIVATE)
target_link_libraries(${MODULE_TARGET} PRIVATE raspakit ${NANOBIND_LIBRARY} ${PYTHON_MODULE_TARGET})
target_include_directories(${MODULE_TARGET} PRIVATE
  ${nanobind_SOURCE_DIR}/include
  ${nanobind_SOURCE_DIR}/ext/robin-map/include
  ${PYTHON_MODULE_INCLUDE_DIRS}
)
include_directories(${Python_INCLUDE_DIRS})
target_compile_options(${MODULE_TARGET} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${MODULE_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(${MODULE_TARGET} PROPERTIES
  PREFIX ""
  SUFFIX "${PYTHON_MODULE_SUFFIX}" # Leverages nanobind's resolved extension macro
)

if(PYPI)

  install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
  if(APPLE)
    install(CODE "
            execute_process(
                COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${MODULE_TARGET}>\"
            )
        ")
  endif()
  install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})

else()

  # the site-packages of the interpreter found above (a bare 'python3' could be another one), on every platform
  if(CONDA AND DEFINED ENV{PYTHON})
    set(RASPALIB_PYTHON_EXECUTABLE "$ENV{PYTHON}")
  elseif(Python_EXECUTABLE)
    set(RASPALIB_PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
  else()
    set(RASPALIB_PYTHON_EXECUTABLE python3)
  endif()

  if(CONDA)
    target_link_options(${MODULE_TARGET} PRIVATE -L$ENV{CONDA_PREFIX}/lib)

    if(DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
      execute_process(
        COMMAND "${RASPALIB_PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
        OUTPUT_VARIABLE RASPALIB_PYTHON_EXT_SUFFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      set_target_properties(${MODULE_TARGET} PROPERTIES SUFFIX "${RASPALIB_PYTHON_EXT_SUFFIX}")
    endif()
  endif()

  execute_process(
    COMMAND "${RASPALIB_PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  install(TARGETS ${MODULE_TARGET} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)

endif()

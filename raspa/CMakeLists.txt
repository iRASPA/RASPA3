project(raspalib)

# The extension module is installed inside the raspalib package and re-exported by its __init__.py, so it
# needs a name of its own: a submodule named raspalib inside a package named raspalib cannot be imported.
set(MODULE_TARGET _raspalib)

if(PYPI)

  if(ABI3)
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
    
    add_library(${MODULE_TARGET} MODULE)
    target_sources(${MODULE_TARGET} PUBLIC raspa_pybindings.cpp)
    target_compile_options(${MODULE_TARGET} PRIVATE -DMODULE_NAME=_raspalib)
    set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
    set_property(TARGET ${MODULE_TARGET} PROPERTY CXX_MODULE_STD ON)
    set_property(TARGET ${MODULE_TARGET} PROPERTY CXX_SCAN_FOR_MODULES ON)
    raspa_target_link_libraries(${MODULE_TARGET} PRIVATE)
    target_link_libraries(${MODULE_TARGET} PRIVATE raspakit nanobind-static-abi3 Python::SABIModule)
    target_include_directories(${MODULE_TARGET} PRIVATE
      ${nanobind_SOURCE_DIR}/include
      ${nanobind_SOURCE_DIR}/ext/robin-map/include
      ${Python_SABI_INCLUDE_DIRS}
    )
    include_directories(${Python_INCLUDE_DIRS})
    target_compile_options(${MODULE_TARGET} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
    set_target_properties(${MODULE_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_target_properties(${MODULE_TARGET} PROPERTIES
      PREFIX ""
      SUFFIX "${NB_SUFFIX_S}" # Leverages nanobind's resolved extension macro
    )
    
    if(WIN32)
      install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
      install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    elseif(LINUX)
      install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
      install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    elseif(APPLE)
        install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
        install(CODE "
                execute_process(
                    COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${MODULE_TARGET}>\"
                )
            ")
        install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
    endif()
  else()
    if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
      find_package(Python REQUIRED COMPONENTS Development)
    else()
      find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
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
    
    nanobind_build_library(nanobind-static)
    
    add_library(${MODULE_TARGET} MODULE)
    target_sources(${MODULE_TARGET} PUBLIC raspa_pybindings.cpp)
    target_compile_options(${MODULE_TARGET} PRIVATE -DMODULE_NAME=_raspalib)
    set_property(SOURCE raspa_pybindings.cpp PROPERTY CXX_SCAN_FOR_MODULES ON)
    set_property(TARGET ${MODULE_TARGET} PROPERTY CXX_MODULE_STD ON)
    set_property(TARGET ${MODULE_TARGET} PROPERTY CXX_SCAN_FOR_MODULES ON)
    raspa_target_link_libraries(${MODULE_TARGET} PRIVATE)
    target_link_libraries(${MODULE_TARGET} PRIVATE raspakit nanobind-static Python::Module)
    target_include_directories(${MODULE_TARGET} PRIVATE
      ${nanobind_SOURCE_DIR}/include
      ${nanobind_SOURCE_DIR}/ext/robin-map/include
      ${Python_INCLUDE_DIRS}
    )
    include_directories(${Python_INCLUDE_DIRS})
    target_compile_options(${MODULE_TARGET} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
    set_target_properties(${MODULE_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_target_properties(${MODULE_TARGET} PROPERTIES
      PREFIX ""
      SUFFIX "${NB_SUFFIX}" # Leverages nanobind's resolved extension macro
    )
    
    if(WIN32)
      install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
      install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    elseif(LINUX)
      install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
      install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
    elseif(APPLE)
        install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/raspalib)
        install(CODE "
                execute_process(
                    COMMAND codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\" \"${SKBUILD_PLATLIB_DIR}/raspalib/$<TARGET_FILE_NAME:${MODULE_TARGET}>\"
                )
            ")
        install(DIRECTORY raspalib DESTINATION ${SKBUILD_PLATLIB_DIR})
    endif()
  endif()

elseif(CONDA)

  if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
    set(Python3_FIND_STRATEGY LOCATION)
    if(NOT Python3_ROOT_DIR)
      set(Python3_ROOT_DIR "$ENV{PREFIX}")
    endif()
    find_package(Python3 REQUIRED COMPONENTS Development)
  else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
  endif()

  if(DEFINED ENV{PYTHON})
    set(RASPALIB_PYTHON_EXECUTABLE "$ENV{PYTHON}")
  else()
    set(RASPALIB_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
  endif()
  
  include(FetchContent)
 
  # explicitly handles backported compatibility patches for GCC 15 system headers
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v3.0.4
  )
  
  FetchContent_GetProperties(pybind11)
  
  if(NOT pybind11_POPULATED)
    FetchContent_MakeAvailable(pybind11)
  endif()

  pybind11_add_module(${MODULE_TARGET} MODULE)
  target_sources(${MODULE_TARGET} PUBLIC raspa_pybindings_pybind11.cpp)
  target_sources(${MODULE_TARGET} PUBLIC FILE_SET CXX_MODULES FILES raspa_pybindings_pybind11.ixx)
  raspa_target_link_libraries(${MODULE_TARGET})
  target_compile_options(${MODULE_TARGET} PRIVATE -DMODULE_NAME=_raspalib)
  target_link_libraries(${MODULE_TARGET} PRIVATE raspakit)
  target_compile_options(${MODULE_TARGET} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
  set_target_properties(${MODULE_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  
  target_link_options(${MODULE_TARGET} PRIVATE -L$ENV{CONDA_PREFIX}/lib)

  execute_process(
    COMMAND "${RASPALIB_PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
    execute_process(
      COMMAND "${RASPALIB_PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
      OUTPUT_VARIABLE RASPALIB_PYTHON_EXT_SUFFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set_target_properties(${MODULE_TARGET} PROPERTIES SUFFIX "${RASPALIB_PYTHON_EXT_SUFFIX}")
  endif()

  install(TARGETS ${MODULE_TARGET} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)

else()

  if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
    find_package(Python3 REQUIRED COMPONENTS Development)
  else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
  endif()
  
  include(FetchContent)
  
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v3.0.4
  )
  
  FetchContent_GetProperties(pybind11)
  
  if(NOT pybind11_POPULATED)
    FetchContent_MakeAvailable(pybind11)
  endif()
  
  pybind11_add_module(${MODULE_TARGET} MODULE)
  target_sources(${MODULE_TARGET} PUBLIC raspa_pybindings_pybind11.cpp)
  target_sources(${MODULE_TARGET} PUBLIC FILE_SET CXX_MODULES FILES raspa_pybindings_pybind11.ixx)
  raspa_target_link_libraries(${MODULE_TARGET})
  target_compile_options(${MODULE_TARGET} PRIVATE -DMODULE_NAME=_raspalib)
  target_link_libraries(${MODULE_TARGET} PRIVATE raspakit)
  target_compile_options(${MODULE_TARGET} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
  set_target_properties(${MODULE_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  
  if(LINUX)
    execute_process(COMMAND python3 -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    install(TARGETS ${MODULE_TARGET} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
    install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
  elseif(WIN32)
  elseif(APPLE)
  endif()

endif()

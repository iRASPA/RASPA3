project(raspalib)

if(UNIX AND DEFINED ENV{CONDA_BUILD_CROSS_COMPILATION})
  find_package(Python3 REQUIRED COMPONENTS Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
endif()

include(FetchContent)

if(CONDA)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v3.0.1
  )
else()
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.10.4
  )
endif()

FetchContent_GetProperties(pybind11)

if(NOT pybind11_POPULATED)
  FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(${PROJECT_NAME} MODULE)
target_sources(${PROJECT_NAME} PUBLIC raspa_pybindings.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES raspa_pybindings.ixx)

raspa_target_link_libraries(${PROJECT_NAME})

target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

if(CONDA)
  target_link_options(${PROJECT_NAME} PRIVATE -L$ENV{CONDA_PREFIX}/lib)
   execute_process(COMMAND python -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS ${PROJECT_NAME} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
elseif(LINUX)
   execute_process(COMMAND python -c "import sysconfig; print(sysconfig.get_path('platlib'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS ${PROJECT_NAME} DESTINATION ${PYTHON_SITE_PACKAGES}/raspalib COMPONENT raspalib)
  install(DIRECTORY raspalib DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT raspalib)
elseif(APPLE)
  install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/lib)
endif()

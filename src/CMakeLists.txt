project(raspakit)

#set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
#set(CMAKE_SHARED_LIBRARY_PREFIX "")

file(GLOB SRC_CPP_FILES 
  ${PROJECT_SOURCE_DIR}/foundationkit/*.cpp
  ${PROJECT_SOURCE_DIR}/mathkit/*.cpp
  ${PROJECT_SOURCE_DIR}/symmetrykit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/properties/*.cpp)

file(GLOB SRC_IXX_FILES 
  ${PROJECT_SOURCE_DIR}/foundationkit/*.ixx
  ${PROJECT_SOURCE_DIR}/mathkit/*.ixx
  ${PROJECT_SOURCE_DIR}/symmetrykit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/properties/*.ixx)


add_library(${PROJECT_NAME} STATIC "")
target_sources(${PROJECT_NAME} PUBLIC ${SRC_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC -stdlib=libc++ $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
add_library(libs::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
#install(TARGETS ${PROJECT_NAME} DESTINATION lib)

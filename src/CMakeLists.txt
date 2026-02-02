project(raspakit)

file(GLOB SRC_CPP_FILES
  ${PROJECT_SOURCE_DIR}/foundationkit/*.cpp
  ${PROJECT_SOURCE_DIR}/mathkit/*.cpp
  ${PROJECT_SOURCE_DIR}/symmetrykit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/rigid_molecules/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/flexible_molecules/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/component/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/system/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/group/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/cross_system/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/intramolecular/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/interpolation/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/properties/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/integrators/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/ruptura/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/grid_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/integration/*.cpp
  ${PROJECT_SOURCE_DIR}/utilitykit/*.cpp)

file(GLOB SRC_IXX_FILES
  ${PROJECT_SOURCE_DIR}/foundationkit/*.ixx
  ${PROJECT_SOURCE_DIR}/mathkit/*.ixx
  ${PROJECT_SOURCE_DIR}/symmetrykit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/rigid_molecules/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/cbmc/flexible_molecules/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/interactions/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/forcefield/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/component/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/system/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/group/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/mc_moves/cross_system/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/intramolecular/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/potentials/interpolation/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/properties/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/integrators/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/ruptura/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/monte_carlo_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/energy_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/grid_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/raspakit/structure_properties/integration/*.ixx
  ${PROJECT_SOURCE_DIR}/utilitykit/*.ixx)

#set(PCH_HEADER "${PROJECT_SOURCE_DIR}/pch.h")
#set(PCH_SOURCE "${PROJECT_SOURCE_DIR}/pch.cpp")
#set(PCH_OUTPUT "${CMAKE_BINARY_DIR}/pch.h.gch")
#
#add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PCH_OUTPUT}
#    COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -x c++-header ${PCH_HEADER} -o ${PCH_OUTPUT}
#    DEPENDS ${PCH_HEADER} ${PCH_SOURCE}
#    COMMENT "Building precompiled header ${PCH_OUTPUT}"
#)

add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PUBLIC ${SRC_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${SRC_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
add_library(libs::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

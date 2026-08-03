project(structurekit)

# The structural analysis of a porous solid: pore diameters and channel dimensionality, void fractions,
# surface areas, pore-size distributions, blocking spheres and tessellations, by geometric, Monte Carlo
# and energy-based routes, on the processor and on the GPU.
#
# This is deliberately a library of its own rather than another directory inside raspakit, and it does not
# depend on raspakit at all. It sits beside it, on the same three kits below both, in the way symmetrykit
# does. What it is given is a cell, a set of atom centres with a type and a charge each, and a table of what
# a pair of types does to one another; a system, a component, a move, a molecule and an interpolation grid
# are all things it never needs to know exist. Making that a link boundary rather than a convention means
# the compiler refuses an engine dependency instead of a reviewer having to notice one.
#
# It is also built once. raspakit is compiled four times over on x86, once per SIMD level, so that the
# Python wheel can dispatch at run time; none of that helps the geometry, graph traversal and OpenCL host
# code here, and none of those variants was ever linked against it.

file(GLOB STRUCTUREKIT_CPP_FILES
  ${PROJECT_SOURCE_DIR}/shared/*.cpp
  ${PROJECT_SOURCE_DIR}/shared_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/diagrams/accessibility/*.cpp
  ${PROJECT_SOURCE_DIR}/diagrams/apollonius/*.cpp
  ${PROJECT_SOURCE_DIR}/diagrams/voronoi/*.cpp
  ${PROJECT_SOURCE_DIR}/diagrams/exact/*.cpp
  ${PROJECT_SOURCE_DIR}/diagrams/grid_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/shared/*.cpp
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/monte_carlo/*.cpp
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/monte_carlo_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/integration/*.cpp
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/integration_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/energy_based/shared/*.cpp
  ${PROJECT_SOURCE_DIR}/energy_based/energy/*.cpp
  ${PROJECT_SOURCE_DIR}/energy_based/energy_opencl/*.cpp
  ${PROJECT_SOURCE_DIR}/brute_force/*.cpp)

file(GLOB STRUCTUREKIT_IXX_FILES
  ${PROJECT_SOURCE_DIR}/shared/*.ixx
  ${PROJECT_SOURCE_DIR}/shared_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/diagrams/accessibility/*.ixx
  ${PROJECT_SOURCE_DIR}/diagrams/apollonius/*.ixx
  ${PROJECT_SOURCE_DIR}/diagrams/voronoi/*.ixx
  ${PROJECT_SOURCE_DIR}/diagrams/exact/*.ixx
  ${PROJECT_SOURCE_DIR}/diagrams/grid_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/shared/*.ixx
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/monte_carlo/*.ixx
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/monte_carlo_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/integration/*.ixx
  ${PROJECT_SOURCE_DIR}/monte_carlo_sampling/integration_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/energy_based/shared/*.ixx
  ${PROJECT_SOURCE_DIR}/energy_based/energy/*.ixx
  ${PROJECT_SOURCE_DIR}/energy_based/energy_opencl/*.ixx
  ${PROJECT_SOURCE_DIR}/brute_force/*.ixx)

add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${STRUCTUREKIT_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${STRUCTUREKIT_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)

# Containers and file formats, vector and matrix arithmetic, and crystal symmetry, and nothing above them.
# The link is PUBLIC because those types appear in the interfaces here, so anything using this library sees
# them too.
target_link_libraries(${PROJECT_NAME} PUBLIC symmetrykit mathkit foundationkit OpenMP::OpenMP_CXX)

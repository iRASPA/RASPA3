project(raspakit)

# The simulation engine: systems, components and molecules, the force field and the interactions between
# atoms, the Monte Carlo moves and the CBMC growth behind them, the dynamics integrators and thermostats,
# minimisation, the equations of state, and the properties and averages a run reports.
#
# It stands on the three kits below it and is what the executables, the Python module and the tests are built
# against. The structural analysis beside it is not part of it and must not become part of it; what it is, is
# something this library is allowed to call. A component asking for its blocking pockets to be worked out from
# the framework rather than listed by hand is the engine asking the analysis a question, and the answer arrives
# before the first molecule is placed, so it has to be asked from in here. The arrow points this way only:
# structurekit still links nothing above the three kits, and the compiler still refuses it if that changes.

file(GLOB RASPAKIT_CPP_FILES
  ${PROJECT_SOURCE_DIR}/*.cpp
  ${PROJECT_SOURCE_DIR}/drivers/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/*.cpp
  ${PROJECT_SOURCE_DIR}/cbmc/flexible_molecules/*.cpp
  ${PROJECT_SOURCE_DIR}/interactions/*.cpp
  ${PROJECT_SOURCE_DIR}/forcefield/*.cpp
  ${PROJECT_SOURCE_DIR}/minimization/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/component/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/system/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/group/*.cpp
  ${PROJECT_SOURCE_DIR}/mc_moves/cross_system/*.cpp
  ${PROJECT_SOURCE_DIR}/potentials/*.cpp
  ${PROJECT_SOURCE_DIR}/potentials/intramolecular/*.cpp
  ${PROJECT_SOURCE_DIR}/potentials/interpolation/*.cpp
  ${PROJECT_SOURCE_DIR}/properties/*.cpp
  ${PROJECT_SOURCE_DIR}/properties_time_evolution/*.cpp
  ${PROJECT_SOURCE_DIR}/integrators/*.cpp
  ${PROJECT_SOURCE_DIR}/utilities/*.cpp)

file(GLOB RASPAKIT_IXX_FILES
  ${PROJECT_SOURCE_DIR}/*.ixx
  ${PROJECT_SOURCE_DIR}/drivers/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/*.ixx
  ${PROJECT_SOURCE_DIR}/cbmc/flexible_molecules/*.ixx
  ${PROJECT_SOURCE_DIR}/interactions/*.ixx
  ${PROJECT_SOURCE_DIR}/forcefield/*.ixx
  ${PROJECT_SOURCE_DIR}/minimization/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/component/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/system/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/group/*.ixx
  ${PROJECT_SOURCE_DIR}/mc_moves/cross_system/*.ixx
  ${PROJECT_SOURCE_DIR}/potentials/*.ixx
  ${PROJECT_SOURCE_DIR}/potentials/intramolecular/*.ixx
  ${PROJECT_SOURCE_DIR}/potentials/interpolation/*.ixx
  ${PROJECT_SOURCE_DIR}/properties/*.ixx
  ${PROJECT_SOURCE_DIR}/properties_time_evolution/*.ixx
  ${PROJECT_SOURCE_DIR}/integrators/*.ixx
  ${PROJECT_SOURCE_DIR}/utilities/*.ixx)

# One library, and no -march anywhere below it: the released packages have to run on any processor of the
# architecture they were built for, so nothing here may assume more than the baseline instruction set.
add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${RASPAKIT_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${RASPAKIT_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME} PUBLIC structurekit symmetrykit mathkit foundationkit)

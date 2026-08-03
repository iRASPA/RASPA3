project(foundationkit)

# The bottom of the stack: containers and index paths, archiving and the HDF5 and JSON wrappers, string and
# character handling, the command-line scanner, the thread pool, the OpenCL context and the unit system.
#
# Nothing here imports anything outside the standard library and this directory, which is why the target
# links against no other kit. Anything that would give it a dependency on vectors, symmetry or the simulation
# engine belongs one layer up rather than here; lammps.ixx was moved out to raspakit for exactly that reason.
file(GLOB FOUNDATIONKIT_CPP_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB FOUNDATIONKIT_IXX_FILES ${PROJECT_SOURCE_DIR}/*.ixx)

add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${FOUNDATIONKIT_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${FOUNDATIONKIT_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)

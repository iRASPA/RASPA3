project(structure_input)

# The adapter between the simulation engine and the structural analysis, and the only thing in the tree that
# is allowed to know both. Anything that wants to measure a framework it has read with the engine's CIF
# reader links this; structurekit itself does not, and cannot.

add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/structure_input.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${PROJECT_SOURCE_DIR}/structure_input.ixx)
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME} PUBLIC raspakit_base structurekit OpenMP::OpenMP_CXX)

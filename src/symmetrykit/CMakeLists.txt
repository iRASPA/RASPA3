project(symmetrykit)

# Crystal symmetry: space groups and their settings, point groups, Seitz matrices and changes of basis, the
# primitive-cell and Niggli reduction, band paths, the CIF and POSCAR parsers, and the Voronoi and Apollonius
# diagrams.
#
# It stands on mathkit for the vector and ring arithmetic and on foundationkit for the parsers' scanner and
# character sets, and knows nothing of the simulation engine above it. structurekit was given the same shape
# deliberately.
file(GLOB SYMMETRYKIT_CPP_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB SYMMETRYKIT_IXX_FILES ${PROJECT_SOURCE_DIR}/*.ixx)

add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${SYMMETRYKIT_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${SYMMETRYKIT_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME} PUBLIC mathkit foundationkit OpenMP::OpenMP_CXX)

# See foundationkit: the SIMD variants of raspakit compile these sources themselves.
set(SYMMETRYKIT_CPP_FILES ${SYMMETRYKIT_CPP_FILES} PARENT_SCOPE)
set(SYMMETRYKIT_IXX_FILES ${SYMMETRYKIT_IXX_FILES} PARENT_SCOPE)

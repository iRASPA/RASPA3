project(mathkit)

# Vectors, matrices and quaternions in two, three and four dimensions, the eigensolvers and the Hessian
# layouts, random numbers, marching cubes, interpolation and the ring arithmetic the symmetry search needs.
#
# It reaches down into foundationkit for archiving and JSON, so that these types can serialise themselves,
# and nowhere else.
file(GLOB MATHKIT_CPP_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB MATHKIT_IXX_FILES ${PROJECT_SOURCE_DIR}/*.ixx)

add_library(${PROJECT_NAME} STATIC ${PCH_HEADER} ${PCH_SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_sources(${PROJECT_NAME} PRIVATE ${MATHKIT_CPP_FILES})
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES ${MATHKIT_IXX_FILES})
target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>)
target_link_libraries(${PROJECT_NAME} PUBLIC foundationkit OpenMP::OpenMP_CXX)

# See foundationkit: the SIMD variants of raspakit compile these sources themselves.
set(MATHKIT_CPP_FILES ${MATHKIT_CPP_FILES} PARENT_SCOPE)
set(MATHKIT_IXX_FILES ${MATHKIT_IXX_FILES} PARENT_SCOPE)

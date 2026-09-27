project(raspa3)

add_executable(${PROJECT_NAME}
        main.cpp
)

target_sources(${PROJECT_NAME} PUBLIC commandline.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES commandline.ixx)

# The commit and compile date shown by --help / --version. The script runs on every build (a git call)
# and rewrites the file only when the commit or the dirty state changes, so it costs a relink only then.
set(BUILD_INFO_CPP ${CMAKE_CURRENT_BINARY_DIR}/build_info.cpp)
add_custom_target(${PROJECT_NAME}_build_info
        COMMAND ${CMAKE_COMMAND}
                -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
                -DOUTPUT=${BUILD_INFO_CPP}
                -DVERSION=${CMAKE_PROJECT_VERSION}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/build_info.cmake
        BYPRODUCTS ${BUILD_INFO_CPP}
        COMMENT "Recording git commit for ${PROJECT_NAME}"
        VERBATIM)
add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_build_info)
target_sources(${PROJECT_NAME} PRIVATE ${BUILD_INFO_CPP})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

raspa_target_link_libraries(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PRIVATE raspakit structurekit)

if(BUILD_MAC_PACKAGE)
        install(CODE "execute_process(COMMAND
    codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
    \"app/${PROJECT_NAME}\")" COMPONENT app)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT app)

project(raspa3-cli)

add_executable(${PROJECT_NAME}
  main.cpp
)

# structure_input is the adapter that turns a framework and a force field into what the structural analysis
# takes. It is compiled in here rather than being a library of its own because this is what needs it, and
# because it is the only place in the tree allowed to know both raspakit and structurekit; the unit tests
# compile their own copy of the same two files.
target_sources(${PROJECT_NAME} PUBLIC commandline.cpp structure_input.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES commandline.ixx structure_input.ixx)

raspa_target_link_libraries(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PRIVATE raspakit structurekit)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
    codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
    \"cli/${PROJECT_NAME}\")" COMPONENT cli)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT cli)

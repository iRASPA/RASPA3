project(raspa3-cli)

add_executable(${PROJECT_NAME}
  main.cpp
)

target_sources(${PROJECT_NAME} PUBLIC commandline.cpp)
target_sources(${PROJECT_NAME} PUBLIC FILE_SET CXX_MODULES FILES commandline.ixx)

raspa_target_link_libraries(${PROJECT_NAME})

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
    codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
    \"cli/${PROJECT_NAME}\")" COMPONENT cli)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT cli)

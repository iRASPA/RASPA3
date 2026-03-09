project(raspa3)

add_executable(${PROJECT_NAME}
        main.cpp
)

raspa_target_link_libraries(${PROJECT_NAME})

if(BUILD_MAC_PACKAGE)
        install(CODE "execute_process(COMMAND
    codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
    \"app/${PROJECT_NAME}\")" COMPONENT app)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT app)

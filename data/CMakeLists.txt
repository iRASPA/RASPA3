install(DIRECTORY forcefields DESTINATION share/raspa3 COMPONENT app)
install(DIRECTORY molecules DESTINATION share/raspa3 COMPONENT app)
install(DIRECTORY manual DESTINATION share/raspa3 COMPONENT manual)

if(WIN32 AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "AMD64")
  install(DIRECTORY dlls/ DESTINATION bin COMPONENT app)
endif()

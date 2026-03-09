project(benchmarks_raspakit)

add_executable(${PROJECT_NAME}
  interactions.cpp
  pow.cpp
  main.cpp)

raspa_target_link_libraries(${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/benchmarks COMPONENT ${PROJECT_NAME})

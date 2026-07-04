project(benchmarks_raspakit)

# interactions.cpp is excluded: it imports the removed 'factory' module and no longer builds
add_executable(${PROJECT_NAME}
  potentials.cpp
  pow.cpp
  main.cpp)

raspa_target_link_libraries(${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PRIVATE
                      raspakit_base
                      benchmark::benchmark)

install(TARGETS ${PROJECT_NAME} DESTINATION share/raspa3/benchmarks COMPONENT ${PROJECT_NAME})

project(raspa)

add_executable(${PROJECT_NAME}
	main.cpp
)

target_link_libraries(${PROJECT_NAME}
	PRIVATE
		libs::raspakit
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

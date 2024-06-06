project(raspa)

add_executable(${PROJECT_NAME}
        main.cpp
)

execute_process(COMMAND bash -c "gfortran --print-file-name libgfortran.a" OUTPUT_VARIABLE LIB_GFORTRAN OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND bash -c "gfortran --print-file-name libquadmath.a" OUTPUT_VARIABLE LIB_QUADMATH OUTPUT_STRIP_TRAILING_WHITESPACE)

target_link_libraries(${PROJECT_NAME}
        PUBLIC
        libs::raspakit
        -static-libstdc++
        -static-libgcc
        OpenMP::OpenMP_CXX
        -l:libc++.a
        -l:libc++abi.a
        -l:liblapack.a
        -l:libblas.a
        -lm
        ${LIB_GFORTRAN}
        ${LIB_QUADMATH}
        -lpthread
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

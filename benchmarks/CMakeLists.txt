include(FetchContent)
FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main) # need main for benchmark::benchmark

if (MSVC)
  set(gtest_force_shared_crt on)     # needed to avoid error LNK2038
endif()

FetchContent_MakeAvailable(googlebenchmark)
  
add_subdirectory(raspakit-benchmarks)

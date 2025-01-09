include(FetchContent)
FetchContent_Declare(
  scran_tests
  GIT_REPOSITORY https://github.com/libscran/scran_tests
  GIT_TAG master
)
FetchContent_MakeAvailable(scran_tests)

add_executable(
    libtest 
    src/aggregate_across_cells.cpp
    src/aggregate_across_genes.cpp
    src/combine_factors.cpp
    src/clean_factor.cpp
)

target_link_libraries(
    libtest
    scran_aggregate
    scran_tests
)

target_compile_options(libtest PRIVATE -Wall -Werror -Wpedantic -Wextra)

option(CODE_COVERAGE "Enable coverage testing" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

option(TEST_DIRTY_INIT "Test dirty vector initialization" OFF)
if(TEST_DIRTY_INIT)
    target_compile_definitions(libtest PRIVATE "SCRAN_AGGREGATE_TEST_INIT=scran_tests::initial_value()")
endif()

include(GoogleTest)
gtest_discover_tests(libtest)

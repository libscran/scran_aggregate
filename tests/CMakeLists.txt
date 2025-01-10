include(FetchContent)
FetchContent_Declare(
  scran_tests
  GIT_REPOSITORY https://github.com/libscran/scran_tests
  GIT_TAG master
)
FetchContent_MakeAvailable(scran_tests)

include(GoogleTest)

option(CODE_COVERAGE "Enable coverage testing" OFF)
set(DO_CODE_COVERAGE OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(DO_CODE_COVERAGE ON)
endif()

macro(decorate_test target)
    target_link_libraries(
        ${target}
        scran_aggregate
        scran_tests
    )

    target_compile_options(${target} PRIVATE -Wall -Werror -Wpedantic -Wextra)

    if(DO_CODE_COVERAGE)
        target_compile_options(${target} PRIVATE -O0 -g --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endif()

    gtest_discover_tests(${target})
endmacro()

add_executable(
    libtest 
    src/aggregate_across_cells.cpp
    src/aggregate_across_genes.cpp
    src/combine_factors.cpp
    src/clean_factor.cpp
)
decorate_test(libtest)

add_executable(
    dirtytest 
    src/aggregate_across_cells.cpp
    src/aggregate_across_genes.cpp
    src/combine_factors.cpp
    src/clean_factor.cpp
)
decorate_test(dirtytest)
target_compile_definitions(dirtytest PRIVATE "SCRAN_AGGREGATE_TEST_INIT=scran_tests::initial_value()")

add_executable(
    cuspartest
    src/aggregate_across_cells.cpp
    src/aggregate_across_genes.cpp
)
decorate_test(cuspartest)
target_compile_definitions(cuspartest PRIVATE SUBPAR_CUSTOM_PARALLELIZE_RANGE=::subpar::test_parallelize_range)

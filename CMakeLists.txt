cmake_minimum_required(VERSION 3.14)

project(aggregate_across_cells
    VERSION 1.0.0
    DESCRIPTION "Aggregate statistics across cells"
    LANGUAGES CXX)

add_library(aggregate_across_cells INTERFACE)
target_compile_features(aggregate_across_cells INTERFACE cxx_std_17)
target_include_directories(aggregate_across_cells INTERFACE include/)

add_subdirectory(extern)
target_link_libraries(aggregate_across_cells INTERFACE tatami::tatami)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
